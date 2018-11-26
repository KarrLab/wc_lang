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
                     Compartment, RateLaw, RateLawExpression, RateLawDirection, SubmodelAlgorithm,
                     BiomassComponent, BiomassReaction, SpeciesTypeType, Parameter,
                     Observable, ObservableExpression, Function, FunctionExpression,
                     ExpressionMethods)
from wc_lang.io import Reader
from wc_lang.prepare import PrepareModel
from wc_lang.expression_utils import ParsedExpression

# configuration
import wc_lang.config.core
config_wc_lang = wc_lang.config.core.get_config()['wc_lang']


class TestPrepareModel(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_model.xlsx')

    def setUp(self):
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)
        self.dfba_submodel = self.model.submodels.get_one(id='submodel_1')
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
