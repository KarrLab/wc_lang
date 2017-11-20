""" Test WC model preparation

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-10-22
:Copyright: 2017, Karr Lab
:License: MIT
"""
import os
import unittest
import six
import numpy as np

from wc_lang.core import (Model, Submodel, ObjectiveFunction, Reaction, SpeciesType, Species,
    Compartment, ReactionParticipant, RateLaw, RateLawEquation, RateLawDirection, SubmodelAlgorithm,
    Concentration, BiomassComponent, BiomassReaction, SpeciesTypeType)
from wc_lang.io import Reader
from wc_lang.prepare import PrepareModel, CheckModel

# configuration
from wc_utils.config.core import ConfigManager
from wc_lang.config import paths as config_paths_wc_lang
config_wc_lang = \
    ConfigManager(config_paths_wc_lang.core).get_config()['wc_lang']


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
                [(1.0, 'reaction_1'), (2.0, 'reaction_2')], [(1.0, 'Metabolism_biomass'),]),
            ('2*Metabolism_biomass + -4.4*reaction_1',
                [(-4.4, 'reaction_1'),], [(2.0, 'Metabolism_biomass'),]),
        ]

        for test_and_result in tests_and_results:

            test, reaction_results, biomass_results = test_and_result
            of.expression = test
            (reactions, biomass_reactions) = parse_dfba_submodel_obj_func(self.dfba_submodel)
            for coeff,reaction in reaction_results:
                self.assertIn((coeff,reaction), reactions)
            self.assertEqual(len(reactions), len(reaction_results))
            for coeff,biomass_reaction in biomass_results:
                self.assertIn((coeff,biomass_reaction), biomass_reactions)
            self.assertEqual(len(biomass_reactions), len(biomass_results))

        error_inputs = [
            ('reaction_1 +', 'Cannot parse'),
            ('reaction_1 * reaction_1', 'Cannot parse'),
            ('reaction_1 + (reaction_1 - reaction_1)', 'Cannot parse'),
            ('reaction_1 + 3*reaction_1', 'Multiple uses'),
            ('reaction_1 , reaction_2', 'Cannot parse'),
            ('x', 'Unknown reaction or biomass reaction id'),
        ]
        for error_input,msg in error_inputs:
            of.expression = error_input
            with self.assertRaises(ValueError) as context:
                parse_dfba_submodel_obj_func(self.dfba_submodel)
            self.assertIn(msg, str(context.exception))

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
            self.prepare_model.apply_default_dfba_submodel_flux_bounds(self.dfba_submodel), (1,1))
        test_non_rev = self.dfba_submodel.reactions.create(
            id='__test_1',
            reversible=False
        )
        test_rev = self.dfba_submodel.reactions.create(
            id='__test_2',
            reversible=True
        )
        self.assertEqual(self.prepare_model.apply_default_dfba_submodel_flux_bounds(self.dfba_submodel),
            (2,2))
        self.prepare_model.apply_default_dfba_submodel_flux_bounds(self.dfba_submodel)
        self.assertEqual(test_non_rev.max_flux, test_rev.max_flux)
        self.assertEqual(-test_rev.min_flux, test_rev.max_flux)
        self.assertEqual(self.prepare_model.apply_default_dfba_submodel_flux_bounds(self.dfba_submodel),
            (0,0))

    def test_run(self):
        of = self.dfba_submodel.objective_function
        of.expression = 'Metabolism_biomass + reaction_1 + reaction_2*2.0'
        self.prepare_model.run()
        self.assertTrue(of.linear)
        of.expression = 'reaction_1*reaction_1'
        self.prepare_model.run()
        self.assertFalse(of.linear)


class TestGapFinding(unittest.TestCase):

    def setUp(self):
        # make model
        self.model = Model(id='model')
        comp = self.model.compartments.create(id='comp')
        self.species = []
        self.num_species = 20
        for i in range(1, self.num_species+1):
            spec_type = self.model.species_types.create(id='spec_type_{}'.format(i),
                type=SpeciesTypeType.metabolite)
            self.species.append(Species(species_type=spec_type, compartment=comp))
        self.dfba_submodel = self.model.submodels.create(
            id='metabolism', algorithm=SubmodelAlgorithm.dfba)

        self.id_idx = 0
        self.prepare_model = PrepareModel(self.model)

    def next_id(self):
        self.id_idx += 1
        return "rxn_{}".format(self.id_idx)

    def make_reaction(self, submodel, reactant, product, reversible=True):
        rxn = submodel.reactions.create(id=self.next_id(), reversible=reversible)
        rxn.participants.create(species=reactant, coefficient=-1)
        rxn.participants.create(species=product, coefficient=1)

    def create_reaction_network(self, submodel, network_type, **kwargs):
        # make networks of reactions with 1 reactant and 1 product
        # first delete all Reactions
        submodel.reactions = []
        if network_type == 'ring':
            # kwargs options: size, reversible
            species = self.species
            if len(species) < kwargs['size']:
                self.fail("not enough species, len(species) < kwargs['size']")
            for r_idx in range(kwargs['size']):
                product_idx = (r_idx+1) % kwargs['size']
                self.make_reaction(submodel, species[r_idx], species[product_idx], kwargs['reversible'])
        else:
            self.Fail("Unknown network type: {}".format(network_type))

    def test_get_inactive_reactions(self):
        # make ring of 3 irreversible reactions
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'size':3, 'reversible':False})

        # no dead end species -> no inactive reactions
        self.assertEqual(self.prepare_model.get_inactive_reactions(self.dfba_submodel, (set(), set())), [])

        # one dead end species -> 2 inactive reactions
        first_specie = self.species[0]
        dead_end_species = set([first_specie])
        inactive_reactions = self.prepare_model.get_inactive_reactions(self.dfba_submodel,
            (set(), dead_end_species))
        self.assertEqual(len(inactive_reactions), 2)
        self.assertIn(self.dfba_submodel.reactions[0], inactive_reactions)
        self.assertIn(self.dfba_submodel.reactions[-1], inactive_reactions)

    def test_find_dead_end_species(self):
        prep_mdl = self.prepare_model

        # make ring of 4 irreversible reactions
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'size':4, 'reversible':False})

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
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'size':4, 'reversible':False})
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
        self.create_reaction_network(self.dfba_submodel, 'ring', **{'size':3, 'reversible':True})
        # ring with first reaction missing -> all species produced and consumed
        del self.dfba_submodel.reactions[0]
        species_not_consumed, species_not_produced = prep_mdl.find_dead_end_species(self.dfba_submodel, set())
        self.assertFalse(species_not_consumed)
        self.assertFalse(species_not_produced)

    def test_identify_dfba_submodel_rxn_gaps(self):
        prep_mdl = self.prepare_model
        size = 4
        kwargs = {'size':size, 'reversible':False}
        # ring of 4 irreversible reactions -> no dead end species or inactive reactions
        self.create_reaction_network(self.dfba_submodel, 'ring', **kwargs)
        (not_consumed, not_produced), inactive_rxns = prep_mdl.identify_dfba_submodel_rxn_gaps(self.dfba_submodel)
        self.assertFalse(not_consumed)
        self.assertFalse(not_produced)
        self.assertFalse(inactive_rxns)

        # ring of 4 irreversible reactions with one missing -> all species dead end and all reactions inactive
        del self.dfba_submodel.reactions[0]
        (not_consumed, not_produced), inactive_rxns = prep_mdl.identify_dfba_submodel_rxn_gaps(self.dfba_submodel)
        species_in_ring = set(self.species[0:size])
        self.assertEqual(not_consumed, species_in_ring)
        self.assertEqual(not_produced, species_in_ring)
        self.assertEqual(sorted(inactive_rxns, key=lambda x: x.id),
            sorted(self.dfba_submodel.reactions, key=lambda x: x.id))

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
        reaction1.min_flux = np.nan
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
        self.assertEquals(len(errors), 1)
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
        rate_law_equation.expression='__ 0'
        self.assertIn("Security risk: rate law expression '__",
            self.check_model.transcode_and_check_rate_law_equations()[0])
        rate_law_equation.expression='not_a_specie[e]'
        self.assertIn("'not_a_specie[e]' not a known specie",
            self.check_model.transcode_and_check_rate_law_equations()[0])
        
        # rate laws that fail evaluation
        rate_law_equation.expression='foo foo'
        self.assertIn("syntax error in transcoded rate law".format(TEST_ID),
            self.check_model.transcode_and_check_rate_law_equations()[0])
        rate_law_equation.expression='cos(0)'
        self.assertIn("name 'cos' is not defined".format(TEST_ID),
            self.check_model.transcode_and_check_rate_law_equations()[0])
        rate_law_equation.expression='{{{*'
        self.assertIn("EOF in multi-line statement",
            self.check_model.transcode_and_check_rate_law_equations()[0])

    def test_verify_reactant_compartments(self):
        for actual,expected in zip(self.check_model.verify_reactant_compartments(), 
            [".*reaction_1 uses specie specie_1 in another compartment: e",
                ".*reaction_1 uses specie specie_2 in another compartment: e",
                ".*'ssa_submodel' must contain a compartment attribute"]):
            six.assertRegex(self, actual, expected)

    def test_run(self):
        self.check_model.run()
