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
    Concentration, BiomassComponent, BiomassReaction)
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
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)
        self.prepare_model = PrepareModel(self.model)

    def test_fill_dfba_submodel_reaction_gaps(self):
        dfba_submodel = Submodel.objects.get_one(id='submodel_1')
        self.assertEqual(
            self.prepare_model.fill_dfba_submodel_reaction_gaps(dfba_submodel), 3)

        produced = set()
        consumed = set()
        GAP_FILLING_RXN_ID_PREFIX = config_wc_lang['GAP_FILLING_RXN_ID_PREFIX']
        GAP_FILLING_RXN_NAME_PREFIX = config_wc_lang['GAP_FILLING_RXN_NAME_PREFIX']
        for rxn in dfba_submodel.reactions:
            if rxn.id.startswith(GAP_FILLING_RXN_ID_PREFIX):
                for part in rxn.participants:
                    if part.coefficient<0:
                        consumed.add(part.species)
                    elif 0<part.coefficient:
                        produced.add(part.species)

        expected_produced = Species.get(['specie_1[e]', 'specie_2[e]'], dfba_submodel.get_species())
        expected_consumed = Species.get(['specie_1[c]'], dfba_submodel.get_species())
        self.assertEqual(produced, set(expected_produced))
        self.assertEqual(consumed, set(expected_consumed))

    def test_confirm_dfba_submodel_obj_func(self):
        dfba_submodel = Submodel.objects.get_one(id='submodel_1')

        confirm_dfba_submodel_obj_func = self.prepare_model.confirm_dfba_submodel_obj_func

        dfba_submodel.algorithm = None
        with self.assertRaises(ValueError) as context:
            confirm_dfba_submodel_obj_func(dfba_submodel)
        self.assertIn("not a dfba submodel", str(context.exception))
        dfba_submodel.algorithm = SubmodelAlgorithm.dfba

        self.assertEqual(confirm_dfba_submodel_obj_func(dfba_submodel), None)

        dfba_submodel.objective_function = None
        self.assertEqual(confirm_dfba_submodel_obj_func(dfba_submodel), None)
        # dfba_submodel should be using its biomass reaction as its objective function
        self.assertEqual(dfba_submodel.objective_function.expression, dfba_submodel.biomass_reaction.id)

        dfba_submodel.objective_function = None
        dfba_submodel.biomass_reaction.id = dfba_submodel.reactions[0].id
        with self.assertRaises(ValueError) as context:
            confirm_dfba_submodel_obj_func(dfba_submodel)

    def test_apply_default_dfba_submodel_flux_bounds(self):
        dfba_submodel = Submodel.objects.get_one(id='submodel_1')
        self.assertEqual(
            self.prepare_model.apply_default_dfba_submodel_flux_bounds(dfba_submodel), (1,1))
        test_non_rev = dfba_submodel.reactions.create(
            id='__test_1',
            reversible=False
        )
        test_rev = dfba_submodel.reactions.create(
            id='__test_2',
            reversible=True
        )
        self.assertEqual(self.prepare_model.apply_default_dfba_submodel_flux_bounds(dfba_submodel),
            (2,2))
        self.prepare_model.apply_default_dfba_submodel_flux_bounds(dfba_submodel)
        self.assertEqual(test_non_rev.max_flux, test_rev.max_flux)
        self.assertEqual(-test_rev.min_flux, test_rev.max_flux)
        self.assertEqual(self.prepare_model.apply_default_dfba_submodel_flux_bounds(dfba_submodel),
            (0,0))

    def test_run(self):
        self.prepare_model.run()


class TestCheckModel(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_check_model_model.xlsx')

    def setUp(self):
        for model in [Submodel, Reaction, SpeciesType]:
            model.objects.reset()
        # read a wc model
        self.model = Reader().run(self.MODEL_FILENAME)
        self.check_model = CheckModel(self.model)

    def test_check_dfba_submodel_1(self):
        dfba_submodel = Submodel.objects.get_one(id='dfba_submodel')
        self.assertEqual(self.check_model.check_dfba_submodel(dfba_submodel), [])

        # delete a reaction's min flux
        reaction1 = Reaction.objects.get_one(id='reaction_1')
        reaction1.min_flux = np.nan
        errors = self.check_model.check_dfba_submodel(dfba_submodel)
        self.assertIn("Error: no min_flux for reaction 'reaction_name_1' in submodel", errors[0])

    def test_check_dfba_submodel_2(self):
        dfba_submodel = Submodel.objects.get_one(id='dfba_submodel')

        # violate reaction.min_flux <= reaction.max_flux
        reaction1 = Reaction.objects.get_one(id='reaction_1')
        reaction1.min_flux = reaction1.max_flux + 1

        # violate reaction.reversible => reaction.min_flux <= 0
        reaction2 = Reaction.objects.get_one(id='reaction_2')
        reaction2.min_flux = 1

        errors = self.check_model.check_dfba_submodel(dfba_submodel)
        self.assertIn("Error: max_flux < min_flux ({} < {}) for reaction '{}' in submodel".format(
            reaction1.max_flux, reaction1.min_flux, reaction1.name), errors[0])
        self.assertIn("Error: 0 < min_flux ({}) for reversible reaction '{}' in submodel".format(
            reaction2.min_flux, reaction2.name), errors[1])

    def test_check_dfba_submodel_3(self):
        dfba_submodel = Submodel.objects.get_one(id='dfba_submodel')

        # remove all BiomassComponents from the BiomassReaction
        dfba_submodel.biomass_reaction.biomass_components = []
        errors = self.check_model.check_dfba_submodel(dfba_submodel)
        self.assertIn("Error: submodel '{}' uses dfba but lacks a biomass reaction".format(dfba_submodel.name),
            errors[0])

        # remove the BiomassReaction
        dfba_submodel.biomass_reaction = None
        errors = self.check_model.check_dfba_submodel(dfba_submodel)
        self.assertIn("Error: submodel '{}' uses dfba but lacks a biomass reaction".format(dfba_submodel.name),
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
