""" Tests of utilities.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang import (Model, Taxon, Submodel, Reaction, SpeciesType, SpeciesTypeType, Species,
                     Compartment, SpeciesCoefficient, BiomassComponent, BiomassReaction,
                     Parameter, Reference, ReferenceType, DatabaseReference, RateLaw,
                     RateLawEquation, SubmodelAlgorithm, Concentration, ObjectiveFunction,
                     Observable, Function, FunctionExpression, StopCondition, StopConditionExpression,
                     ObservableExpression)
from wc_lang import util
import shutil
import tempfile
import unittest


class TestUtil(unittest.TestCase):
    """ Test utilities """

    def setUp(self):
        self.model = mdl = Model()

        self.comp_0 = comp_0 = mdl.compartments.create(id='comp_0', name='compartment 0')
        self.comp_1 = comp_1 = mdl.compartments.create(id='comp_1', name='compartment 1')
        self.compartments = compartments = [comp_0, comp_1]

        self.species_types = species_types = []
        self.species = species = []
        self.concentrations = concentrations = []
        for i in range(8):
            spec_type = mdl.species_types.create(id='spec_type_{}'.format(
                i), name='species type {}'.format(i), type=SpeciesTypeType.metabolite)
            species_types.append(spec_type)

            if i != 3:
                spec = Species(species_type=spec_type, compartment=comp_0)
            else:
                spec = Species(species_type=spec_type, compartment=comp_1)
            species.append(spec)

            conc = Concentration(species=spec, value=1)
            concentrations.append(conc)

        self.submdl_0 = submdl_0 = mdl.submodels.create(id='submdl_0', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_1 = submdl_1 = mdl.submodels.create(id='submdl_1', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_2 = submdl_2 = mdl.submodels.create(id='submdl_2', algorithm=SubmodelAlgorithm.dfba)
        self.submodels = [submdl_0, submdl_1, submdl_2]

        self.rxn_0 = rxn_0 = submdl_0.reactions.create(id='rxn_0')
        rxn_0.participants.create(species=species[0], coefficient=-2)
        rxn_0.participants.create(species=species[1], coefficient=-3)
        rxn_0.participants.create(species=species[2], coefficient=1)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[5].get_primary_attribute()),
            modifiers=species[5:6])
        rate_law_0 = rxn_0.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_1 = rxn_1 = submdl_1.reactions.create(id='rxn_1')
        rxn_1.participants.create(species=species[0], coefficient=-2)
        rxn_1.participants.create(species=species[1], coefficient=-3)
        rxn_1.participants.create(species=species[3], coefficient=2)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[6].get_primary_attribute()),
            modifiers=species[6:7])
        rate_law_1 = rxn_1.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_2 = rxn_2 = submdl_2.reactions.create(id='rxn_2')
        rxn_2.participants.create(species=species[0], coefficient=-2)
        rxn_2.participants.create(species=species[1], coefficient=-3)
        rxn_2.participants.create(species=species[4], coefficient=1)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[7].get_primary_attribute()),
            modifiers=species[7:8])
        rate_law_2 = rxn_2.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.reactions = [rxn_0, rxn_1, rxn_2]
        self.rate_laws = [rate_law_0, rate_law_1, rate_law_2]

        self.parameters = parameters = []
        self.references = references = []
        self.database_references = database_references = []
        for i in range(3):
            param = mdl.parameters.create(id='param_{}'.format(i))
            parameters.append(param)

            ref = param.references.create(id='ref_{}'.format(i), type=ReferenceType.misc)
            references.append(ref)

            x_ref = ref.database_references.create(database='Y', id='x')
            database_references.append(x_ref)

    def test_get_model_size(self):
        model = self.model
        size = util.get_model_size(model)
        self.assertEqual(3, size['submodels'])
        self.assertEqual(8, size['species_types'])
        self.assertEqual(8, size['species'])
        self.assertEqual(3, size['reactions'])
        self.assertEqual(2, size['compartments'])
        self.assertEqual(3, size['parameters'])
        self.assertEqual(3, size['references'])

    def test_get_model_summary(self):
        model = self.model
        summary = util.get_model_summary(model)
        self.assertIsInstance(summary, str)

    def test_get_reaction_string(self):
        species_types = self.species_types
        species = self.species

        self.assertIn(util.get_reaction_string(self.rxn_0), [
            '[{0}]: ({1}) {2} + ({3}) {4} ==> {5}'.format(self.comp_0.id, 2,
                                                          species_types[0].id, 3, species_types[1].id, species_types[2].id),
            '[{0}]: ({3}) {4} + ({1}) {2} ==> {5}'.format(self.comp_0.id, 2,
                                                          species_types[0].id, 3, species_types[1].id, species_types[2].id),
        ])

        self.assertIn(util.get_reaction_string(self.rxn_1), [
            '({0}) {1} + ({2}) {3} ==> (2) {4}'.format(2,
                                                       species[0].serialize(), 3, species[1].serialize(), species[3].serialize()),
            '({2}) {3} + ({0}) {1} ==> (2) {4}'.format(2,
                                                       species[0].serialize(), 3, species[1].serialize(), species[3].serialize()),
        ])

    def test_get_models(self):
        non_inline_models = set([
            Model, Taxon,
            Submodel, Compartment, SpeciesType, Observable, Concentration,
            Reaction, RateLaw, BiomassComponent, BiomassReaction, Parameter,
            Function, StopCondition, Reference, DatabaseReference,
        ])
        inline_models = set([
            Species, SpeciesCoefficient, RateLawEquation, ObjectiveFunction,
            FunctionExpression, StopConditionExpression, ObservableExpression
        ])
        self.assertEqual(set(util.get_models()), non_inline_models | inline_models)
        self.assertEqual(set(util.get_models(inline=False)), non_inline_models)

    def test_set_git_repo_metadata_from_path(self):
        model = Model()
        self.assertEqual(model.url, '')

        util.set_git_repo_metadata_from_path(model, path='.')
        self.assertIn(model.url, [
            'https://github.com/KarrLab/wc_lang.git',
            'ssh://git@github.com/KarrLab/wc_lang.git',
            'git@github.com:KarrLab/wc_lang.git',
        ])

    def test_set_git_repo_metadata_from_path_error(self):
        tempdir = tempfile.mkdtemp()

        model = Model()
        self.assertEqual(model.url, '')

        with self.assertRaisesRegexp(ValueError, 'is not a Git repository'):
            util.set_git_repo_metadata_from_path(model, path=tempdir)
        self.assertEqual(model.url, '')

        shutil.rmtree(tempdir)
