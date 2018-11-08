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

from wc_lang import (Model, Submodel, ObjectiveFunction, Reaction, SpeciesType, Species,
                     Compartment, RateLaw, RateLawEquation, RateLawDirection, SubmodelAlgorithm,
                     BiomassComponent, BiomassReaction, SpeciesTypeType, Parameter,
                     Observable, ObservableExpression, Function, FunctionExpression,
                     ExpressionMethods)
from wc_lang.io import Reader
from wc_lang.prepare import PrepareModel, CheckModel, AnalyzeModel
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
        with self.assertRaisesRegex(ValueError, ' not a dfba submodel$'):
            self.prepare_model.create_dfba_exchange_rxns(Submodel(algorithm=SubmodelAlgorithm.ssa), None)

    def test_confirm_dfba_submodel_obj_func(self):

        confirm_dfba_submodel_obj_func = self.prepare_model.confirm_dfba_submodel_obj_func

        self.assertEqual(confirm_dfba_submodel_obj_func(self.dfba_submodel), None)

        self.dfba_submodel.algorithm = None
        with self.assertRaises(ValueError) as context:
            confirm_dfba_submodel_obj_func(self.dfba_submodel)
        self.assertIn("not a dfba submodel", str(context.exception))
        self.dfba_submodel.algorithm = SubmodelAlgorithm.dfba

        self.dfba_submodel.objective_function = None
        self.assertEqual(confirm_dfba_submodel_obj_func(self.dfba_submodel), None)
        # self.dfba_submodel should be using its biomass reaction as its objective function
        self.assertEqual(self.dfba_submodel.objective_function.expression,
                         self.dfba_submodel.biomass_reaction.id)
        self.assertEqual(self.dfba_submodel.objective_function.reactions, [])
        self.assertEqual(self.dfba_submodel.objective_function.biomass_reaction_coefficients, [1.0])

    def test_parse_dfba_submodel_obj_func(self):
        parse_dfba_submodel_obj_func = self.prepare_model.parse_dfba_submodel_obj_func

        of = self.dfba_submodel.objective_function

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
            of.expression = test
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
            of.expression = error_input
            with self.assertRaises(ValueError) as context:
                parse_dfba_submodel_obj_func(self.dfba_submodel)
            self.assertIn(msg, str(context.exception))

        # test exception
        with self.assertRaisesRegex(ValueError, ' not a dfba submodel$'):
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
        of = self.dfba_submodel.objective_function
        of.expression = 'Metabolism_biomass + reaction_1 + reaction_2*2.0'
        (reactions, biomass_reactions) = self.prepare_model.parse_dfba_submodel_obj_func(self.dfba_submodel)
        PrepareModel.assign_linear_objective_fn(self.dfba_submodel, reactions, biomass_reactions)
        self.assertEqual(of.biomass_reactions[0].id, 'Metabolism_biomass')
        self.assertEqual(of.biomass_reaction_coefficients[0], 1.0)
        coeffs_n_ids = zip(of.reaction_coefficients, [r.id for r in of.reactions])
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
        with self.assertRaisesRegex(ValueError, ' not a dfba submodel$'):
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
        of = self.dfba_submodel.objective_function
        of.expression = 'Metabolism_biomass + reaction_1 + reaction_2*2.0'
        self.prepare_model.run()
        self.assertTrue(of.linear)
        of.expression = 'reaction_1*reaction_1'
        self.prepare_model.run()
        self.assertFalse(of.linear)


class TestAnalyzeModel(unittest.TestCase):

    RNX_ID_PREFIX = 'rxn'
    SPECIES_ID_PREFIX = 'spec_type'
    default_max_flux = 10000

    def rxn_id(self, n):
        return "{}_{}".format(TestAnalyzeModel.RNX_ID_PREFIX, n)

    def sp_id(self, n):
        return "{}_{}".format(TestAnalyzeModel.SPECIES_ID_PREFIX, n)

    def next_id(self):
        self.id_idx += 1
        return self.rxn_id(self.id_idx)

    def setUp(self):
        # make model
        self.model = Model(id='model')
        comp = self.model.compartments.create(id='comp')
        self.species = []
        self.num_species = 20
        for i in range(1, self.num_species+1):
            spec_type = self.model.species_types.create(id=self.sp_id(i),
                                                        type=SpeciesTypeType.metabolite)
            self.species.append(Species(species_type=spec_type, compartment=comp))
        self.dfba_submodel = self.model.submodels.create(
            id='metabolism', algorithm=SubmodelAlgorithm.dfba)

        self.id_idx = 0
        self.analyze_model = AnalyzeModel(self.model)

    def make_reaction(self, submodel, reactant, product, **kwargs):
        reversible = True
        if 'reversible' in kwargs:
            reversible = kwargs['reversible']
        max_flux = TestAnalyzeModel.default_max_flux
        if 'max_flux' in kwargs:
            max_flux = kwargs['max_flux']
        rxn = submodel.reactions.create(id=self.next_id(), reversible=reversible, max_flux=max_flux)
        rxn.participants.create(species=reactant, coefficient=-1)
        rxn.participants.create(species=product, coefficient=1)

    def create_reaction_network(self, submodel, network_type, **kwargs):
        # make networks of reactions with 1 reactant and 1 product
        # first delete all Reactions
        submodel.reactions = []
        if network_type == 'ring':
            # kwargs options: num_rxn, reversible, max_flux
            species = self.species
            if len(species) < kwargs['num_rxn']:
                self.fail("not enough species, len(species) < kwargs['num_rxn']")
            for reactant_idx in range(kwargs['num_rxn']):
                product_idx = (reactant_idx+1) % kwargs['num_rxn']
                self.make_reaction(submodel, species[reactant_idx], species[product_idx], **kwargs)
        else:
            self.Fail("Unknown network type: {}".format(network_type))

    def test_get_inactive_reactions(self):
        # make ring of 3 irreversible reactions
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'num_rxn': 3, 'reversible': False})

        # no dead end species -> no inactive reactions
        self.assertEqual(self.analyze_model.get_inactive_reactions(self.dfba_submodel, (set(), set())), [])

        # one dead end species -> 2 inactive reactions
        first_specie = self.species[0]
        dead_end_species = set([first_specie])
        inactive_reactions = self.analyze_model.get_inactive_reactions(self.dfba_submodel,
                                                                       (set(), dead_end_species))
        self.assertEqual(len(inactive_reactions), 2)
        self.assertIn(self.dfba_submodel.reactions[0], inactive_reactions)
        self.assertIn(self.dfba_submodel.reactions[-1], inactive_reactions)

    def test_find_dead_end_species(self):
        prep_mdl = self.analyze_model

        # make ring of 4 irreversible reactions
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'num_rxn': 4, 'reversible': False})

        # ring with no inactive reactions -> no dead end species
        species_not_consumed, species_not_produced = prep_mdl.find_dead_end_species(self.dfba_submodel, set())
        self.assertFalse(species_not_consumed)
        self.assertFalse(species_not_produced)

        # ring with first reaction missing ->
        #   species_not_consumed = first reaction's reactant
        #   species_not_produced = first reaction's product
        for part in self.dfba_submodel.reactions[0].participants:
            if part.coefficient == -1:
                reactant = part.species
            if part.coefficient == 1:
                product = part.species
        del self.dfba_submodel.reactions[0]
        species_not_consumed, species_not_produced = prep_mdl.find_dead_end_species(self.dfba_submodel, set())
        self.assertEqual(species_not_consumed.pop(), reactant)
        self.assertEqual(species_not_produced.pop(), product)
        self.assertFalse(species_not_consumed)
        self.assertFalse(species_not_produced)

        # make ring of 4 irreversible reactions
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'num_rxn': 4, 'reversible': False})
        # ring with first reaction inactive ->
        #   species_not_consumed = first reaction's reactant
        #   species_not_produced = first reaction's product
        species_not_consumed, species_not_produced = prep_mdl.find_dead_end_species(self.dfba_submodel,
                                                                                    set([self.dfba_submodel.reactions[0]]))
        self.assertEqual(species_not_consumed.pop(), reactant)
        self.assertEqual(species_not_produced.pop(), product)
        self.assertFalse(species_not_consumed)
        self.assertFalse(species_not_produced)

        # make ring of reversible reactions
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'num_rxn': 3, 'reversible': True})
        # ring with first reaction missing -> all species produced and consumed
        del self.dfba_submodel.reactions[0]
        species_not_consumed, species_not_produced = prep_mdl.find_dead_end_species(self.dfba_submodel, set())
        self.assertFalse(species_not_consumed)
        self.assertFalse(species_not_produced)

    def test_identify_dfba_submodel_rxn_gaps(self):
        prep_mdl = self.analyze_model
        num_rxn = 4
        kwargs = {'num_rxn': num_rxn, 'reversible': False}
        # ring of 4 irreversible reactions -> no dead end species or inactive reactions
        self.create_reaction_network(self.dfba_submodel, 'ring', **kwargs)
        (not_consumed, not_produced), inactive_rxns = prep_mdl.identify_dfba_submodel_rxn_gaps(self.dfba_submodel)
        self.assertFalse(not_consumed)
        self.assertFalse(not_produced)
        self.assertFalse(inactive_rxns)

        # ring of 4 irreversible reactions with one missing -> all species dead end and all reactions inactive
        del self.dfba_submodel.reactions[0]
        (not_consumed, not_produced), inactive_rxns = prep_mdl.identify_dfba_submodel_rxn_gaps(self.dfba_submodel)
        species_in_ring = set(self.species[0:num_rxn])
        self.assertEqual(not_consumed, species_in_ring)
        self.assertEqual(not_produced, species_in_ring)
        self.assertEqual(sorted(inactive_rxns, key=lambda x: x.id),
                         sorted(self.dfba_submodel.reactions, key=lambda x: x.id))

        # check exceptions
        with self.assertRaisesRegex(ValueError, 'not a dfba submodel'):
            prep_mdl.identify_dfba_submodel_rxn_gaps(Submodel(algorithm=SubmodelAlgorithm.ssa))

    def test_digraph_of_rxn_network(self):
        self.run_test_on_digraph_of_rxn_network(5, False)

    def test_digraph_of_rxn_network_reversible(self):
        self.run_test_on_digraph_of_rxn_network(5, True)

    def run_test_on_digraph_of_rxn_network(self, num_rxn, reversible):
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'num_rxn': num_rxn, 'reversible': reversible})
        g = self.analyze_model.digraph_of_rxn_network(self.dfba_submodel)

        for n in g.nodes():
            if isinstance(n, Reaction):
                self.assertTrue(TestAnalyzeModel.RNX_ID_PREFIX in n.id)
            elif isinstance(n, Species):
                self.assertTrue(TestAnalyzeModel.SPECIES_ID_PREFIX in n.id())

        # test expected vs. actual edges
        # expected edge id pairs
        expected_edges = set()
        # forward:
        for i in range(1, num_rxn+1):
            rxn_2_sp_edge = (self.rxn_id(i), self.sp_id((i % num_rxn)+1))
            expected_edges.add(rxn_2_sp_edge)
            sp_2_rxn_edge = (self.sp_id(i), self.rxn_id(i))
            expected_edges.add(sp_2_rxn_edge)
        if reversible:
            for i in range(1, num_rxn+1):
                rxn_2_sp_edge = (self.rxn_id(i), self.sp_id(i))
                expected_edges.add(rxn_2_sp_edge)
                sp_2_rxn_edge = (self.sp_id((i % num_rxn)+1), self.rxn_id(i))
                expected_edges.add(sp_2_rxn_edge)

        graph_edges = set()
        for s, d in g.edges():
            ids = []
            for n in [s, d]:
                if isinstance(n, Reaction):
                    ids.append(n.id)
                if isinstance(n, Species):
                    # remove compartment suffix '[some_comp]'
                    sp_type_id = n.id().split('[')[0]
                    ids.append(sp_type_id)
            s_id, d_id = ids
            graph_edges.add((s_id, d_id))
        self.assertEqual(expected_edges, graph_edges)

    def test_unbounded_paths(self):
        num_rxn = 8

        # irrreversible reactions
        # unbounded network
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'num_rxn': num_rxn, 'reversible': False,
                                                                    'max_flux': float('inf')})
        path_len = 2*num_rxn-1
        g = self.analyze_model.digraph_of_rxn_network(self.dfba_submodel)
        paths = self.analyze_model.unbounded_paths(g, self.species[0], [self.species[num_rxn-1]])
        self.assertEqual(len(paths), 1)
        self.assertEqual(len(paths[0]), path_len)

        # bounded network
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'num_rxn': num_rxn, 'reversible': False})
        g = self.analyze_model.digraph_of_rxn_network(self.dfba_submodel)
        paths = self.analyze_model.unbounded_paths(g, self.species[0],
                                                   [self.species[num_rxn-1]], min_non_finite_ub=self.default_max_flux+1)
        self.assertEqual(len(paths), 0)

        # reversible reactions, paths on both sides of ring
        # unbounded network
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'num_rxn': num_rxn, 'reversible': True,
                                                                    'max_flux': float('inf')})
        g = self.analyze_model.digraph_of_rxn_network(self.dfba_submodel)
        paths = self.analyze_model.unbounded_paths(g, self.species[0], [self.species[num_rxn//2]])
        self.assertEqual(len(paths), 2)
        for p in paths:
            self.assertEqual(len(p), num_rxn+1)

        # bounded network
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'num_rxn': num_rxn, 'reversible': True})
        g = self.analyze_model.digraph_of_rxn_network(self.dfba_submodel)
        paths = self.analyze_model.unbounded_paths(g, self.species[0], [self.species[num_rxn//2]],
                                                   min_non_finite_ub=self.default_max_flux+1)
        self.assertEqual(len(paths), 0)

        # test exceptions
        with self.assertRaisesRegex(ValueError, "'ex_species' should be a Species instance, but "):
            self.analyze_model.unbounded_paths(None, 'species', None)

        with self.assertRaisesRegex(ValueError, "elements of 'obj_fn_species' should be Species instances, but "):
            self.analyze_model.unbounded_paths(None, Species(), ['species'])

    def test_path_bounds_analysis(self):
        # read a wc model
        MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_model.xlsx')
        Submodel.objects.reset()
        Reaction.objects.reset()
        BiomassReaction.objects.reset()
        self.model = Reader().run(MODEL_FILENAME)
        self.dfba_submodel = Submodel.objects.get_one(id='submodel_1')
        self.analyze_model = AnalyzeModel(self.model)

        for rxn in self.dfba_submodel.reactions:
            rxn.max_flux = 0
        paths = self.analyze_model.path_bounds_analysis(self.dfba_submodel)
        for k in paths.keys():
            self.assertEqual(paths[k], [])

        for rxn in self.dfba_submodel.reactions:
            rxn.max_flux = float('inf')
        paths = self.analyze_model.path_bounds_analysis(self.dfba_submodel)
        self.assertEqual(len(paths['specie_1[e]']), 2)


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
        self.dfba_submodel.biomass_reaction.biomass_components = []
        errors = self.check_model.check_dfba_submodel(self.dfba_submodel)
        self.assertIn("Error: submodel '{}' uses dfba but lacks a biomass reaction".format(self.dfba_submodel.name),
                      errors[0])

        # remove the BiomassReaction
        self.dfba_submodel.biomass_reaction = None
        errors = self.check_model.check_dfba_submodel(self.dfba_submodel)
        self.assertIn("Error: submodel '{}' uses dfba but lacks a biomass reaction".format(self.dfba_submodel.name),
                      errors[0])

        # remove the objective function
        self.dfba_submodel.objective_function = None
        errors = self.check_model.check_dfba_submodel(self.dfba_submodel)
        self.assertIn("Error: submodel '{}' uses dfba but lacks an objective function".format(self.dfba_submodel.name),
                      errors[0])

    def test_check_dfba_submodel_4(self):

        # remove a reaction to test that all species used in biomass reactions are defined
        del self.dfba_submodel.reactions[-1]
        errors = self.check_model.check_dfba_submodel(self.dfba_submodel)
        self.assertEqual(len(errors), 1)
        six.assertRegex(self, errors[0],
                        "Error: undefined species '.*' in biomass reaction '.*' used by submodel")

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

    def test_transcode_and_check_rate_law_equations(self):
        # good laws
        self.assertEqual(self.check_model.transcode_and_check_rate_law_equations(), [])

        # test errors
        # redefine one reaction
        rate_law_equation = RateLawEquation(
            expression='',
            transcoded='',
        )
        rate_law = RateLaw(
            equation=rate_law_equation,
        )
        rate_law_equation.rate_law = rate_law
        a_reaction = self.model.get_reactions().pop()
        a_reaction.rate_laws = [rate_law]
        TEST_ID = 'test_id'
        a_reaction.id = TEST_ID

        # rate laws that fail transcoding
        rate_law_equation.expression = '__ 0'
        self.assertIn("Security risk: rate law expression '__",
                      self.check_model.transcode_and_check_rate_law_equations()[0])
        rate_law_equation.expression = 'not_a_specie[e]'
        self.assertIn("'not_a_specie[e]' not a known specie",
                      self.check_model.transcode_and_check_rate_law_equations()[0])
        rate_law_equation.expression = 'foo foo'
        self.assertIn("'foo' not a known parameter".format(TEST_ID),
                      self.check_model.transcode_and_check_rate_law_equations()[0])
        rate_law_equation.expression = 'cos(0)'
        self.assertIn("'cos' not a known parameter".format(TEST_ID),
                      self.check_model.transcode_and_check_rate_law_equations()[0])

        # rate laws that fail evaluation                
        rate_law_equation.expression = '{{{*'
        self.assertIn("EOF in multi-line statement",
                      self.check_model.transcode_and_check_rate_law_equations()[0])

        submodel = self.model.submodels.get_one(algorithm=SubmodelAlgorithm.ssa)
        submodel.reactions[1].rate_laws[0].equation = None
        self.check_model.transcode_and_check_rate_law_equations()

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
            Observable:{}
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
            obs_expr, e = ExpressionMethods.make_expression_obj(Observable, 'obs_{}'.format((i+1)%num), objects)
            self.assertEqual(e, None)
            id = "obs_{}".format(i)
            objects[Observable][id].expression = obs_expr
        errors = self.check_model.verify_acyclic_dependencies([Observable])
        self.assertEqual(len(errors), 1)
        for e in errors:
            self.assertIn('dependency cycle among Observables', e)

        # Function
        objects = {
            Parameter:{},
            Observable:{},
            Function:{}
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
