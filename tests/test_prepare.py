""" Test WC model preparation

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-10-22
:Copyright: 2017, Karr Lab
:License: MIT
"""
import ast
import mock
import numpy
import os
import re
import six
import unittest

from wc_lang import (Model, Submodel, Reaction, SpeciesType, Species,
                     Compartment, RateLaw, RateLawEquation, RateLawDirection, SubmodelAlgorithm,
                     BiomassComponent, BiomassReaction, SpeciesTypeType, Parameter,
                     Observable, ObservableExpression, Function, FunctionExpression,
                     ExpressionMethods)
from wc_lang.io import Reader
from wc_lang.prepare import PrepareModel, CheckModel
from wc_lang.expression_utils import WcLangExpression

# configuration
import wc_lang.config.core
config_wc_lang = wc_lang.config.core.get_config()['wc_lang']


class TestPrepareModel(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_model.xlsx')

    def setUp(self):
        Submodel.objects.reset()
        Reaction.objects.reset()
        BiomassReaction.objects.reset()
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)
        self.dfba_submodel = Submodel.objects.get_one(id='submodel_1')
        self.prepare_model = PrepareModel(self.model)
        self.id_idx = 0

    def test_create_dfba_exchange_rxns(self):
        EXTRACELLULAR_COMPARTMENT_ID = config_wc_lang['EXTRACELLULAR_COMPARTMENT_ID']

        self.assertEqual(
            self.prepare_model.create_dfba_exchange_rxns(self.dfba_submodel, EXTRACELLULAR_COMPARTMENT_ID), 2)

        # should add these exchange reactions:
        # -> specie_1[e]
        # -> specie_2[e]
        EXCHANGE_RXN_ID_PREFIX = config_wc_lang['EXCHANGE_RXN_ID_PREFIX']
        species_found = set()
        for rxn in self.dfba_submodel.reactions:
            if EXCHANGE_RXN_ID_PREFIX in rxn.id:
                self.assertEqual(-float('inf'), rxn.min_flux)
                self.assertEqual(float('inf'), rxn.max_flux)
                self.assertEqual(1, len(rxn.participants))
                for participant in rxn.participants:
                    self.assertEqual(1, participant.coefficient)
                    species_found.add(participant.species)

        self.assertEqual(species_found,
                         set(Species.get(['specie_1[e]', 'specie_2[e]'], self.dfba_submodel.get_species())))

        # test exception
        with self.assertRaisesRegex(ValueError, ' not a dFBA submodel$'):
            self.prepare_model.create_dfba_exchange_rxns(Submodel(algorithm=SubmodelAlgorithm.ssa), None)

    def test_confirm_dfba_submodel_obj_func(self):

        confirm_dfba_submodel_obj_func = self.prepare_model.confirm_dfba_submodel_obj_func

        self.assertEqual(confirm_dfba_submodel_obj_func(self.dfba_submodel), None)

        self.dfba_submodel.algorithm = None
        with self.assertRaises(ValueError) as context:
            confirm_dfba_submodel_obj_func(self.dfba_submodel)
        self.assertIn("not a dFBA submodel", str(context.exception))
        self.dfba_submodel.algorithm = SubmodelAlgorithm.dfba

        self.dfba_submodel.dfba_obj = None
        self.assertEqual(confirm_dfba_submodel_obj_func(self.dfba_submodel), None)
        # self.dfba_submodel should be using its biomass reaction as its objective function
        self.assertEqual(self.dfba_submodel.dfba_obj.expression.expression,
                         self.dfba_submodel.biomass_reactions[0].id)
        self.assertEqual(self.dfba_submodel.dfba_obj.expression.reactions, [])
        self.assertEqual(self.dfba_submodel.dfba_obj.expression.biomass_reaction_coefficients, [1.0])

    def test_parse_dfba_submodel_obj_func(self):
        parse_dfba_submodel_obj_func = self.prepare_model.parse_dfba_submodel_obj_func

        of = self.dfba_submodel.dfba_obj

        # list of tests and expected results
        # (test, [(coef, reaction_id), ... ] [(coeff, biomass_reaction_id), ... ])
        tests_and_results = [
            ('reaction_1',
                [(1.0, 'reaction_1')], []),
            ('Metabolism_biomass + (reaction_1 + reaction_2*2.0)',
                [(1.0, 'reaction_1'), (2.0, 'reaction_2')], [(1.0, 'Metabolism_biomass'), ]),
            ('2*Metabolism_biomass + -4.4*reaction_1',
                [(-4.4, 'reaction_1'), ], [(2.0, 'Metabolism_biomass'), ]),
        ]

        for test_and_result in tests_and_results:

            test, reaction_results, biomass_results = test_and_result
            of.expression.expression = test
            (reactions, biomass_reactions) = parse_dfba_submodel_obj_func(self.dfba_submodel)
            for coeff, reaction in reaction_results:
                self.assertIn((coeff, reaction), reactions)
            self.assertEqual(len(reactions), len(reaction_results))
            for coeff, biomass_reaction in biomass_results:
                self.assertIn((coeff, biomass_reaction), biomass_reactions)
            self.assertEqual(len(biomass_reactions), len(biomass_results))

        error_inputs = [
            ('reaction_1 +', 'Cannot parse'),
            ('reaction_1 * reaction_1', 'Cannot parse'),
            ('reaction_1 + (reaction_1 - reaction_1)', 'Cannot parse'),
            ('reaction_1 + 3*reaction_1', 'Multiple uses'),
            ('reaction_1 , reaction_2', 'Cannot parse'),
            ('x', 'Unknown reaction or biomass reaction id'),
        ]
        for error_input, msg in error_inputs:
            of.expression.expression = error_input
            with self.assertRaises(ValueError) as context:
                parse_dfba_submodel_obj_func(self.dfba_submodel)
            self.assertIn(msg, str(context.exception))

        # test exception
        with self.assertRaisesRegex(ValueError, ' not a dFBA submodel$'):
            self.prepare_model.parse_dfba_submodel_obj_func(Submodel(algorithm=SubmodelAlgorithm.ssa))

    def test__proc_mult(self):
        node = ast.UnaryOp()
        node.left = ast.UnaryOp()
        node.right = ast.UnaryOp()
        node.left.op = ast.USub()
        node.right.op = ast.USub()
        node.left.operand = ast.Num()
        node.right.operand = ast.Num()
        with self.assertRaisesRegex(ValueError, "bad Mult"):
            self.prepare_model._proc_mult(node, [])

    def test_assign_linear_objective_fn(self):
        of = self.dfba_submodel.dfba_obj
        of.expression.expression = 'Metabolism_biomass + reaction_1 + reaction_2*2.0'
        (reactions, biomass_reactions) = self.prepare_model.parse_dfba_submodel_obj_func(self.dfba_submodel)
        PrepareModel.assign_linear_objective_fn(self.dfba_submodel, reactions, biomass_reactions)
        self.assertEqual(of.expression.biomass_reactions[0].id, 'Metabolism_biomass')
        self.assertEqual(of.expression.biomass_reaction_coefficients[0], 1.0)
        coeffs_n_ids = zip(of.expression.reaction_coefficients,
                           [r.id for r in of.expression.reactions])
        self.assertIn((1.0, 'reaction_1'), coeffs_n_ids)
        self.assertIn((2.0, 'reaction_2'), coeffs_n_ids)

    def test_apply_default_dfba_submodel_flux_bounds(self):
        self.assertEqual(
            self.prepare_model.apply_default_dfba_submodel_flux_bounds(self.dfba_submodel), (1, 1))
        test_non_rev = self.dfba_submodel.reactions.create(
            id='__test_1',
            reversible=False
        )
        test_rev = self.dfba_submodel.reactions.create(
            id='__test_2',
            reversible=True
        )
        self.assertEqual(self.prepare_model.apply_default_dfba_submodel_flux_bounds(self.dfba_submodel),
                         (2, 2))
        self.prepare_model.apply_default_dfba_submodel_flux_bounds(self.dfba_submodel)
        self.assertEqual(test_non_rev.max_flux, test_rev.max_flux)
        self.assertEqual(-test_rev.min_flux, test_rev.max_flux)
        self.assertEqual(self.prepare_model.apply_default_dfba_submodel_flux_bounds(self.dfba_submodel),
                         (0, 0))

        # test exception
        with self.assertRaisesRegex(ValueError, ' not a dFBA submodel$'):
            self.prepare_model.apply_default_dfba_submodel_flux_bounds(Submodel(algorithm=SubmodelAlgorithm.ssa))

        submodel = Submodel(algorithm=SubmodelAlgorithm.dfba)
        submodel.reactions.create(min_flux=float('nan'))
        with mock.patch('wc_lang.prepare.config_wc_lang', {}):
            with self.assertRaisesRegex(ValueError, "cannot obtain default_min_flux_bound and default_max_flux_bound="):
                self.prepare_model.apply_default_dfba_submodel_flux_bounds(submodel)

        submodel = Submodel(algorithm=SubmodelAlgorithm.dfba)
        submodel.reactions.create(min_flux=float('nan'))
        with mock.patch('wc_lang.prepare.config_wc_lang', {'default_min_flux_bound': 1, 'default_max_flux_bound': -1}):
            with self.assertRaisesRegex(ValueError, "default flux bounds violate 0 <= default_min_flux_bound <= default_max_flux_bound:"):
                self.prepare_model.apply_default_dfba_submodel_flux_bounds(submodel)

    def test_run(self):
        of = self.dfba_submodel.dfba_obj

        of.expression.expression = 'Metabolism_biomass + reaction_1 + reaction_2*2.0'
        self.prepare_model.run()
        self.assertEqual(set([rxn.id for rxn in of.expression.reactions]),
                         set(['reaction_1', 'reaction_2']))
        self.assertEqual(set([rxn.id for rxn in of.expression.biomass_reactions]),
                         set(['Metabolism_biomass']))
        self.assertTrue(of.expression.linear)

        of.expression.expression = 'reaction_1*reaction_1'
        self.prepare_model.run()
        self.assertFalse(of.expression.linear)


class TestCheckModel(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_check_model_model.xlsx')

    def setUp(self):
        for model in [Submodel, Reaction, SpeciesType]:
            model.objects.reset()
        # read a wc model
        self.model = Reader().run(self.MODEL_FILENAME)
        self.dfba_submodel = Submodel.objects.get_one(id='dfba_submodel')
        self.check_model = CheckModel(self.model)

    def test_check_dfba_submodel_1(self):
        self.assertEqual(self.check_model.check_dfba_submodel(self.dfba_submodel), [])

        # delete a reaction's min flux
        reaction1 = Reaction.objects.get_one(id='reaction_1')
        reaction1.min_flux = numpy.nan
        errors = self.check_model.check_dfba_submodel(self.dfba_submodel)
        self.assertIn("Error: no min_flux for reaction 'reaction_name_1' in submodel", errors[0])

    def test_check_dfba_submodel_2(self):

        # violate reaction.min_flux <= reaction.max_flux
        reaction1 = Reaction.objects.get_one(id='reaction_1')
        reaction1.min_flux = reaction1.max_flux + 1

        # violate reaction.reversible => reaction.min_flux <= 0
        reaction2 = Reaction.objects.get_one(id='reaction_2')
        reaction2.min_flux = 1

        errors = self.check_model.check_dfba_submodel(self.dfba_submodel)
        self.assertIn("Error: max_flux < min_flux ({} < {}) for reaction '{}' in submodel".format(
            reaction1.max_flux, reaction1.min_flux, reaction1.name), errors[0])
        self.assertIn("Error: 0 < min_flux ({}) for reversible reaction '{}' in submodel".format(
            reaction2.min_flux, reaction2.name), errors[1])

    def test_check_dfba_submodel_3(self):

        # remove all BiomassComponents from the BiomassReaction
        self.dfba_submodel.biomass_reactions[0].biomass_components = []
        errors = self.check_model.check_dfba_submodel(self.dfba_submodel)
        self.assertIn("Error: submodel '{}' uses dFBA but lacks a biomass reaction".format(self.dfba_submodel.name),
                      errors[0])

        # remove the BiomassReaction
        self.dfba_submodel.biomass_reaction = None
        errors = self.check_model.check_dfba_submodel(self.dfba_submodel)
        self.assertIn("Error: submodel '{}' uses dFBA but lacks a biomass reaction".format(self.dfba_submodel.name),
                      errors[0])

        # remove the objective function
        self.dfba_submodel.dfba_obj = None
        errors = self.check_model.check_dfba_submodel(self.dfba_submodel)
        self.assertIn("Error: submodel '{}' uses dFBA but lacks an objective function".format(self.dfba_submodel.name),
                      errors[0])

    def test_check_dynamic_submodel(self):
        ssa_submodel = Submodel.objects.get_one(id='ssa_submodel')
        self.assertEqual(self.check_model.check_dynamic_submodel(ssa_submodel), [])

        reaction_4 = Reaction.objects.get_one(id='reaction_4')
        # add reaction_4 backward ratelaw -> not reversible but has backward error
        reaction_4_ratelaw = reaction_4.rate_laws[0]
        reaction_4.rate_laws[0].direction = RateLawDirection.backward
        errors = self.check_model.check_dynamic_submodel(ssa_submodel)
        self.assertIn("is not reversible but has a 'backward' rate law specified", errors[0])

        # remove reaction_4 forward ratelaw -> no rate law error
        reaction_4.rate_laws = []
        errors = self.check_model.check_dynamic_submodel(ssa_submodel)
        self.assertIn("has no rate law specified", errors[0])

        # put back the good rate law for reaction_4
        reaction_4_ratelaw.direction = RateLawDirection.forward
        reaction_4.rate_laws = [reaction_4_ratelaw]
        self.assertEqual(self.check_model.check_dynamic_submodel(ssa_submodel), [])

        # remove reaction_3 backward ratelaw -> reversible but only forward error
        reaction_3 = Reaction.objects.get_one(id='reaction_3')
        del reaction_3.rate_laws[1:]
        errors = self.check_model.check_dynamic_submodel(ssa_submodel)
        self.assertIn("is reversible but has only a 'forward' rate law specified", errors[0])

    def test_run(self):
        self.check_model.run()

    def test_run_exception(self):
        self.dfba_submodel.reactions[0].min_flux = float('nan')
        with self.assertRaisesRegex(ValueError, 'no min_flux'):
            CheckModel(self.model).run()

    def test_verify_species_types(self):
        self.assertEqual(self.check_model.verify_species_types(), [])
        SpeciesType.objects.get_one(id='specie_4').molecular_weight = float('NaN')
        SpeciesType.objects.get_one(id='specie_6').molecular_weight = -1.
        expected_errors = [
            "species types must contain positive molecular weights, but the MW for specie_4 is nan",
            "species types must contain positive molecular weights, but the MW for specie_6 is -1.0"
        ]
        actual_errors = self.check_model.verify_species_types()
        self.assertEqual(set(expected_errors), set(actual_errors))

    def make_some_objects(self):
        objects = {}
        objects[Parameter] = {}
        objects[Observable] = {}
        objects[Function] = {}
        for i in range(1, 4):
            id = "fun_{}".format(i)
            objects[Function][id] = self.model.functions.create(id=id)
            id = "obs_{}".format(i)
            objects[Observable][id] = self.model.observables.create(id=id)
        return objects

    def test_verify_acyclic_dependencies(self):
        # Observable
        objects = {
            Observable: {}
        }
        num = 3
        for i in range(num):
            id = "obs_{}".format(i)
            objects[Observable][id] = ExpressionMethods.make_obj(self.model, Observable, id, '', {},
                                                                 allow_invalid_objects=True)
        errors = self.check_model.verify_acyclic_dependencies([Observable])
        self.assertEqual(errors, [])

        # create cyclic dependencies: each observable references the next, mod num
        for i in range(num):
            obs_expr, e = ExpressionMethods.make_expression_obj(Observable, 'obs_{}'.format((i+1) % num), objects)
            self.assertEqual(e, None)
            id = "obs_{}".format(i)
            objects[Observable][id].expression = obs_expr
        errors = self.check_model.verify_acyclic_dependencies([Observable])
        self.assertEqual(len(errors), 1)
        for e in errors:
            self.assertIn('dependency cycle among Observables', e)

        # Function
        objects = {
            Parameter: {},
            Observable: {},
            Function: {}
        }

        id = 'fun_3'
        # First 2 functions call the next one; make before referencing
        objects[Function][id] = ExpressionMethods.make_obj(self.model, Function, id, '', {},
                                                           allow_invalid_objects=True)
        for i in range(2, 0, -1):
            id = 'fun_{}'.format(i)
            objects[Function][id] = \
                ExpressionMethods.make_obj(self.model, Function, id, 'fun_{}()'.format(i+1), objects)
        errors = self.check_model.verify_acyclic_dependencies([Function])
        self.assertEqual(errors, [])

        # test cyclic recursion: have fun_1 call itself
        fun_expr, error = ExpressionMethods.make_expression_obj(Function, 'fun_1() + fun_2()', objects)
        self.assertEqual(error, None)
        objects[Function]['fun_1'].expression = fun_expr
        # and fun_3 call fun_1
        objects[Function]['fun_3'].expression = \
            ExpressionMethods.make_expression_obj(Function, 'fun_1()', objects)[0]

        errors = self.check_model.verify_acyclic_dependencies([Function])
        self.assertEqual(len(errors), 2)
        for e in errors:
            self.assertIn('dependency cycle among Functions', e)
