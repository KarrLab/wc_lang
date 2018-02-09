""" Tests of input/output.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang.core import (Model, Taxon, TaxonRank, Submodel, ObjectiveFunction, Reaction, SpeciesType, SpeciesTypeType,
                          Species, Compartment, ReactionParticipant, BiomassComponent, BiomassReaction,
                          Parameter, Reference, ReferenceType, DatabaseReference,
                          RateLaw, RateLawEquation, SubmodelAlgorithm, Concentration)
from wc_lang.io import Writer, Reader, convert, create_template
from wc_utils.workbook.io import read as read_workbook
import obj_model.io
import os
import shutil
import tempfile
import unittest


class TestCreateTemplate(unittest.TestCase):

    def setUp(self):
        _, self.filename = tempfile.mkstemp(suffix='.xlsx')

    def tearDown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)

    def test_create_template(self):
        create_template(self.filename)
        self.assertEqual(Reader().run(self.filename), None)


class TestSimpleModel(unittest.TestCase):

    def setUp(self):
        self.model = mdl = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1b')

        mdl.taxon = Taxon(id='taxon', name='test taxon', rank=TaxonRank.species)

        self.comp_0 = comp_0 = mdl.compartments.create(id='comp_0', name='compartment 0', initial_volume=1.25)
        self.comp_1 = comp_1 = mdl.compartments.create(id='comp_1', name='compartment 1', initial_volume=2.5)
        self.compartments = compartments = [comp_0, comp_1]

        self.species_types = species_types = []
        self.species = species = []
        self.concentrations = concentrations = []
        for i in range(8):
            spec_type = mdl.species_types.create(
                id='spec_type_{}'.format(i),
                name='species type {}'.format(i),
                type=SpeciesTypeType.metabolite,
                structure='C' * (i + 1),
                empirical_formula='C' + str(i + 1),
                molecular_weight=12 * (i + 1),
                charge=i + 1)
            species_types.append(spec_type)

            if i != 3:
                spec = Species(species_type=spec_type, compartment=comp_0)
            else:
                spec = Species(species_type=spec_type, compartment=comp_1)
            species.append(spec)

            conc = Concentration(species=spec, value=3 * i)
            concentrations.append(conc)

        self.submdl_0 = submdl_0 = mdl.submodels.create(
            id='submodel_0', name='submodel 0', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_1 = submdl_1 = mdl.submodels.create(
            id='submodel_1', name='submodel 1', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_2 = submdl_2 = mdl.submodels.create(
            id='submodel_2', name='submodel 2', algorithm=SubmodelAlgorithm.dfba)
        self.submodels = submodels = [submdl_0, submdl_1, submdl_2]

        participants = {}
        def get_or_create_participant(species=None, coefficient=None):
            part_serialized = ReactionParticipant._serialize(species, coefficient)
            if part_serialized not in participants:
                participants[part_serialized] = ReactionParticipant(species=species, coefficient=coefficient)
            return participants[part_serialized]

        self.rxn_0 = rxn_0 = submdl_0.reactions.create(id='rxn_0', name='reaction 0')

        rxn_0.participants.append(get_or_create_participant(species=species[0], coefficient=-2))
        rxn_0.participants.append(get_or_create_participant(species=species[1], coefficient=-3))
        rxn_0.participants.append(get_or_create_participant(species=species[2], coefficient=1))
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[5].serialize()),
            modifiers=species[5:6])
        rate_law_0 = rxn_0.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_1 = rxn_1 = submdl_1.reactions.create(id='rxn_1', name='reaction 1')
        rxn_1.participants.append(get_or_create_participant(species=species[0], coefficient=-2))
        rxn_1.participants.append(get_or_create_participant(species=species[1], coefficient=-3))
        rxn_1.participants.append(get_or_create_participant(species=species[3], coefficient=2))
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[6].serialize()),
            modifiers=species[6:7])
        rate_law_1 = rxn_1.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_2 = rxn_2 = submdl_2.reactions.create(id='rxn_2', name='reaction 2')
        rxn_2.participants.append(get_or_create_participant(species=species[0], coefficient=-2))
        rxn_2.participants.append(get_or_create_participant(species=species[1], coefficient=-3))
        rxn_2.participants.append(get_or_create_participant(species=species[4], coefficient=1))
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[7].serialize()),
            modifiers=species[7:8])
        rate_law_2 = rxn_2.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.reactions = [rxn_0, rxn_1, rxn_2]
        self.rate_laws = [rate_law_0, rate_law_1, rate_law_2]

        self.parameters = parameters = []
        self.references = references = []
        self.database_references = database_references = []
        for i in range(3):
            param = mdl.parameters.create(
                id='param_{}'.format(i), name='parameter {}'.format(i),
                value=i * 4, units='dimensionless')
            param.submodels = submodels[i:i + 1]
            parameters.append(param)

            ref = param.references.create(
                id='ref_{}'.format(i), name='reference {}'.format(i),
                type=ReferenceType.misc)
            references.append(ref)

            x_ref = ref.database_references.create(database='x', id='y' * (i + 1),
                                                url='http://x.com/{}'.format('y' * (i + 1)))
            database_references.append(x_ref)

        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_write_read(self):
        filename = os.path.join(self.dirname, 'model.xlsx')

        Writer().run(filename, self.model)
        model = Reader().run(filename)
        self.assertEqual(model.validate(), None)

        self.assertTrue(model.is_equal(self.model))
        self.assertEqual(self.model.difference(model), '')

    def test_convert(self):
        filename_xls1 = os.path.join(self.dirname, 'model1.xlsx')
        filename_xls2 = os.path.join(self.dirname, 'model2.xlsx')
        filename_csv = os.path.join(self.dirname, 'model-*.csv')

        Writer().run(filename_xls1, self.model)

        convert(filename_xls1, filename_csv)
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'model-Model.csv')))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'model-Taxon.csv')))
        model = Reader().run(filename_csv)
        self.assertTrue(model.is_equal(self.model))

        convert(filename_csv, filename_xls2)
        model = Reader().run(filename_xls2)
        self.assertTrue(model.is_equal(self.model))


class TestExampleModel(unittest.TestCase):

    def setUp(self):
        _, self.filename = tempfile.mkstemp(suffix='.xlsx')

    def tearDown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)

    def test_read_write(self):
        fixture_filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'example-model.xlsx')

        model = Reader().run(fixture_filename)
        self.assertEqual(model.validate(), None)

        # compare excel files
        Writer().run(self.filename, model)
        original = read_workbook(fixture_filename)
        copy = read_workbook(self.filename)
        # note that models must be sorted by id for this assertion to hold
        self.assertEqual(copy, original)

        # compare models
        model2 = Reader().run(self.filename)
        self.assertTrue(model2.is_equal(model))
        self.assertTrue(model.difference(model2) == '')


class TestReaderException(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test(self):
        model1 = Model(id='model1', name='test model', version='0.0.1a', wc_lang_version='0.0.1')
        model2 = Model(id='model2', name='test model', version='0.0.1a', wc_lang_version='0.0.1')
        filename = os.path.join(self.tempdir, 'model.xlsx')
        obj_model.io.Writer().run(filename, [model1, model2], Writer.model_order)

        with self.assertRaisesRegexp(ValueError, ' should only define one model$'):
            Reader().run(filename)
