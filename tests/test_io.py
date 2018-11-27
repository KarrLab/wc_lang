""" Tests of input/output.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2016-11-10
:Copyright: 2016-2018, Karr Lab
:License: MIT
"""

from wc_lang import (Model, Taxon, TaxonRank, Submodel, Reaction, SpeciesType, SpeciesTypeType,
                     Species, Compartment, SpeciesCoefficient,
                     DfbaNetComponent, DfbaNetReaction,
                     Parameter, Reference, ReferenceType, DatabaseReference, Function, FunctionExpression,
                     StopCondition, StopConditionExpression,
                     Observable, ObservableExpression,
                     RateLaw, RateLawExpression, RateLawDirection,
                     SubmodelAlgorithm, Concentration, ConcentrationUnit,
                     DfbaObjective, DfbaObjectiveExpression)
from wc_lang import io
from wc_lang.io import Writer, Reader, convert, create_template
from wc_utils.workbook.io import read as read_workbook, write as write_workbook
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
        create_template(self.filename, set_repo_metadata_from_path=False)
        self.assertIsInstance(Reader().run(self.filename), Model)


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
            spec.id = Species.gen_id(spec.species_type.id, spec.compartment.id)
            spec.model = mdl
            species.append(spec)

            conc = Concentration(id=Concentration.gen_id(spec.id),
                                 model=mdl,
                                 species=spec, value=3 * i, units=ConcentrationUnit.M)
            concentrations.append(conc)

        species_coefficients = {}

        def get_or_create_species_coefficient(species=None, coefficient=None):
            part_serialized = SpeciesCoefficient._serialize(species, coefficient)
            if part_serialized not in species_coefficients:
                species_coefficients[part_serialized] = SpeciesCoefficient(species=species, coefficient=coefficient)
            return species_coefficients[part_serialized]

        self.observables = observables = []
        objects = {Species: {}, Observable: {}}
        for i in range(8):
            expr_parts = []
            for j in range(i + 1):
                objects[Species][species[j].id] = species[j]
                expr_parts.append("{} * {}".format(j + 1, species[j].id))
            obs_expr, _ = ObservableExpression.deserialize(' + '.join(expr_parts), objects)
            id = 'obs_{}'.format(i)
            obs = mdl.observables.create(id=id, expression=obs_expr)
            observables.append(obs)
            objects[Observable][id] = obs
        for i in range(3):
            expr_parts = []
            for j in range(i + 1):
                expr_parts.append("{} * {}".format(j + 1, observables[j].id))
            obs_expr, _ = ObservableExpression.deserialize(' + '.join(expr_parts), objects)
            obs = mdl.observables.create(id='obs_{}'.format(i + 8), expression=obs_expr)
            observables.append(obs)

        objects = {Observable: {o.id: o for o in observables}}
        self.functions = functions = []
        for i in range(8):
            obs_expr = ' + 2 * '.join(o.id for o in observables[0:i + 1])
            func_expr, _ = FunctionExpression.deserialize(obs_expr,
                                                          objects)
            func = mdl.functions.create(id='func_{}'.format(i),
                                        expression=func_expr)
            functions.append(func)

        self.submdl_0 = submdl_0 = mdl.submodels.create(
            id='submodel_0', name='submodel 0', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_1 = submdl_1 = mdl.submodels.create(
            id='submodel_1', name='submodel 1', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_2 = submdl_2 = mdl.submodels.create(
            id='submodel_2', name='submodel 2', algorithm=SubmodelAlgorithm.dfba)
        self.submodels = submodels = [submdl_0, submdl_1, submdl_2]

        self.rxn_0 = rxn_0 = submdl_0.reactions.create(
            id='rxn_0', name='reaction 0', model=mdl)

        rxn_0.participants.append(get_or_create_species_coefficient(species=species[0], coefficient=-2))
        rxn_0.participants.append(get_or_create_species_coefficient(species=species[1], coefficient=-3))
        rxn_0.participants.append(get_or_create_species_coefficient(species=species[2], coefficient=1))
        expression = RateLawExpression(
            expression='k_cat_0 * {0} / (k_m_0 + {0})'.format(species[5].serialize()),
            modifiers=species[5:6])
        expression.parameters.create(id='k_cat_0', value=2, model=mdl)
        expression.parameters.create(id='k_m_0', value=1, model=mdl)
        rate_law_0 = rxn_0.rate_laws.create(
            id=RateLaw.gen_id(rxn_0.id, RateLawDirection.forward.name),
            model=mdl,
            direction=RateLawDirection.forward,
            expression=expression)

        self.rxn_1 = rxn_1 = submdl_1.reactions.create(
            id='rxn_1', name='reaction 1', model=mdl)
        rxn_1.participants.append(get_or_create_species_coefficient(species=species[0], coefficient=-2))
        rxn_1.participants.append(get_or_create_species_coefficient(species=species[1], coefficient=-3))
        rxn_1.participants.append(get_or_create_species_coefficient(species=species[3], coefficient=2))
        expression = RateLawExpression(
            expression='k_cat_1 * {0} / (k_m_1 + {0})'.format(species[6].serialize()),
            modifiers=species[6:7])
        expression.parameters.create(id='k_cat_1', value=2, model=mdl)
        expression.parameters.create(id='k_m_1', value=1, model=mdl)
        rate_law_1 = rxn_1.rate_laws.create(
            id=RateLaw.gen_id(rxn_1.id, RateLawDirection.forward.name),
            model=mdl,
            direction=RateLawDirection.forward,
            expression=expression)

        self.rxn_2 = rxn_2 = submdl_2.reactions.create(
            id='rxn_2', name='reaction 2', model=mdl)
        rxn_2.participants.append(get_or_create_species_coefficient(species=species[0], coefficient=-2))
        rxn_2.participants.append(get_or_create_species_coefficient(species=species[1], coefficient=-3))
        rxn_2.participants.append(get_or_create_species_coefficient(species=species[4], coefficient=1))
        expression = RateLawExpression(
            expression='k_cat_2 * {0} / (k_m_2 + {0})'.format(species[7].serialize()),
            modifiers=species[7:8])
        expression.parameters.create(id='k_cat_2', value=2, model=mdl)
        expression.parameters.create(id='k_m_2', value=1, model=mdl)
        rate_law_2 = rxn_2.rate_laws.create(
            id=RateLaw.gen_id(rxn_2.id, RateLawDirection.forward.name),
            model=mdl,
            direction=RateLawDirection.forward,
            expression=expression)

        submdl_2.dfba_obj = DfbaObjective(id=DfbaObjective.gen_id('submodel_2'), model=mdl)
        submdl_2.dfba_obj.expression = DfbaObjectiveExpression(expression='rxn_2', reactions=[rxn_2])

        self.reactions = [rxn_0, rxn_1, rxn_2]
        self.rate_laws = [rate_law_0, rate_law_1, rate_law_2]

        self.parameters = parameters = []
        self.references = references = []
        self.db_refs = db_refs = []
        for i in range(3):
            param = mdl.parameters.create(
                id='param_{}'.format(i), name='parameter {}'.format(i),
                value=i * 4, units='dimensionless')
            param.submodels = submodels[i:i + 1]
            parameters.append(param)

            ref = param.references.create(
                id='ref_{}'.format(i), name='reference {}'.format(i),
                type=ReferenceType.misc)
            ref.model = mdl
            references.append(ref)

            x_ref = ref.db_refs.create(database='x', id='y' * (i + 1))
            db_refs.append(x_ref)

        self.stop_conditions = stop_conditions = []
        for i in range(3):
            cond = mdl.stop_conditions.create(id='stop_cond_{}'.format(i))
            expr = '{} > 1'.format(' + '.join(o.id for o in observables[0:i+1]))
            cond.expression = StopConditionExpression.deserialize(expr, objects)[0]
            self.stop_conditions.append(cond)

        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_write_read(self):
        filename = os.path.join(self.dirname, 'model.xlsx')

        Writer().run(self.model, filename, set_repo_metadata_from_path=False)
        model = Reader().run(filename)
        self.assertEqual(model.validate(), None)

        self.assertTrue(model.is_equal(self.model))
        self.assertEqual(self.model.difference(model), '')

    def test_write_with_repo_md(self):
        _, filename = tempfile.mkstemp(suffix='.xlsx', dir='.')

        self.assertEqual(self.model.url, '')

        Writer().run(self.model, filename, set_repo_metadata_from_path=True)
        self.assertIn(self.model.url, [
            'https://github.com/KarrLab/wc_lang.git',
            'ssh://git@github.com/KarrLab/wc_lang.git',
            'git@github.com:KarrLab/wc_lang.git',
        ])

        os.remove(filename)

    def test_write_read_sloppy(self):
        filename = os.path.join(self.dirname, 'model.xlsx')

        Writer().run(self.model, filename, set_repo_metadata_from_path=False)

        wb = read_workbook(filename)
        row = wb['Model'].pop(0)
        wb['Model'].insert(1, row)
        write_workbook(filename, wb)

        with self.assertRaisesRegex(ValueError, "The columns of worksheet 'Model' must be defined in this order"):
            Reader().run(filename)
        model = Reader().run(filename, strict=False)
        self.assertEqual(model.validate(), None)

        self.assertTrue(model.is_equal(self.model))
        self.assertEqual(self.model.difference(model), '')

    def test_convert(self):
        filename_xls1 = os.path.join(self.dirname, 'model1.xlsx')
        filename_xls2 = os.path.join(self.dirname, 'model2.xlsx')
        filename_csv = os.path.join(self.dirname, 'model-*.csv')

        Writer().run(self.model, filename_xls1, set_repo_metadata_from_path=False)

        convert(filename_xls1, filename_csv)
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'model-Model.csv')))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'model-Taxon.csv')))
        model = Reader().run(filename_csv)
        self.assertTrue(model.is_equal(self.model))

        convert(filename_csv, filename_xls2)
        model = Reader().run(filename_xls2)
        self.assertTrue(model.is_equal(self.model))

    def test_convert_sloppy(self):
        filename_xls1 = os.path.join(self.dirname, 'model1.xlsx')
        filename_xls2 = os.path.join(self.dirname, 'model2.xlsx')
        filename_csv = os.path.join(self.dirname, 'model-*.csv')

        Writer().run(self.model, filename_xls1, set_repo_metadata_from_path=False)

        wb = read_workbook(filename_xls1)
        row = wb['Model'].pop(0)
        wb['Model'].insert(1, row)
        write_workbook(filename_xls1, wb)

        with self.assertRaisesRegex(ValueError, "The columns of worksheet 'Model' must be defined in this order"):
            convert(filename_xls1, filename_csv)
        convert(filename_xls1, filename_csv, strict=False)

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
        Writer().run(model, self.filename, set_repo_metadata_from_path=False)
        original = read_workbook(fixture_filename)
        copy = read_workbook(self.filename)
        # note that models must be sorted by id for this assertion to hold
        for sheet in original.keys():
            for i_row, (copy_row, original_row) in enumerate(zip(copy[sheet], original[sheet])):
                self.assertEqual(copy_row, original_row,
                                 msg='Rows {} of {} sheets are not equal'.format(i_row, sheet))
            self.assertEqual(copy[sheet], original[sheet], msg='{} sheets are not equal'.format(sheet))

        self.assertEqual(copy, original)

        # compare models
        model2 = Reader().run(self.filename)
        self.assertTrue(model2.is_equal(model))
        self.assertTrue(model.difference(model2) == '')

    def test_rate_law_expressions_with_multiple_model_types(self):
        fixture_filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'example-model.xlsx')
        model = Reader().run(fixture_filename)
        rate_laws = model.get_rate_laws(id='AK_AMP-backward')
        self.assertEqual(len(rate_laws), 1)
        rate_law = rate_laws[0]
        self.assertEqual(len(rate_law.expression.modifiers), 1)
        self.assertEqual(rate_law.expression.modifiers[0].id, 'Adk_Protein[c]')
        self.assertEqual(len(rate_law.expression.parameters), 1)
        self.assertEqual(rate_law.expression.parameters[0].id, 'k_cat_ak')
        self.assertEqual(len(rate_law.expression.observables), 1)
        self.assertEqual(rate_law.expression.observables[0].id, 'AXP_c')
        self.assertEqual(len(rate_law.expression.functions), 1)
        self.assertEqual(rate_law.expression.functions[0].id, 'func_1')


class TestReaderException(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test(self):
        model1 = Model(id='model1', name='test model', version='0.0.1a', wc_lang_version='0.0.1')
        model2 = Model(id='model2', name='test model', version='0.0.1a', wc_lang_version='0.0.1')
        filename = os.path.join(self.tempdir, 'model.xlsx')
        obj_model.io.WorkbookWriter().run(filename, [model1, model2], Writer.model_order, include_all_attributes=False)

        with self.assertRaisesRegex(ValueError, ' should define one model$'):
            Reader().run(filename)


class TestReadNoModel(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test(self):
        filename = os.path.join(self.tempdir, 'model.xlsx')
        obj_model.io.WorkbookWriter().run(filename, [], io.Writer.model_order, include_all_attributes=False)
        self.assertEqual(Reader().run(filename), None)


class ImplicitRelationshipsTestCase(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_write_parameter(self):
        model = Model(id='model', version='0.0.1', wc_lang_version='0.0.1')
        submodel = model.submodels.create(id='submodel')
        reaction = submodel.reactions.create(id='reaction', model=model)
        rate_law = reaction.rate_laws.create(model=model, direction=RateLawDirection.forward)
        rate_law_eq = rate_law.expression = RateLawExpression(expression='parameter')
        parameter = rate_law_eq.parameters.create(id='parameter', model=model)

        filename = os.path.join(self.tempdir, 'model.xlsx')
        Writer().run(model, filename, set_repo_metadata_from_path=False)

        parameter.model = Model(id='model2', version='0.0.1', wc_lang_version='0.0.1')
        with self.assertRaisesRegex(ValueError, 'must be set to the instance of `Model`'):
            Writer().run(model, filename, set_repo_metadata_from_path=False)

    def test_write_other(self):
        model = Model(id='model', version='0.0.1', wc_lang_version='0.0.1')
        species_type = model.species_types.create(id='species_type')
        compartment = model.compartments.create(id='compartment')
        species = Species(
            id=Species.gen_id(species_type.id, compartment.id),
            model=model,
            species_type=species_type,
            compartment=compartment)
        s_id = species.serialize()
        obs_expr, _ = ObservableExpression.deserialize(s_id, {Species: {s_id: species}})
        observable = model.observables.create(id='observable', expression=obs_expr)

        filename = os.path.join(self.tempdir, 'model.xlsx')
        Writer().run(model, filename, set_repo_metadata_from_path=False)

        model2 = Model(id='model2', version='0.0.1', wc_lang_version='0.0.1')
        observable.model = model2
        with self.assertRaisesRegex(ValueError, 'must be set to the instance of `Model`'):
            Writer().run(model, filename, set_repo_metadata_from_path=False)

    def test_read(self):
        filename = os.path.join(self.tempdir, 'model.xlsx')
        obj_model.io.WorkbookWriter().run(filename, [Submodel(id='submodel')], Writer.model_order, include_all_attributes=False)
        with self.assertRaisesRegex(ValueError, 'cannot contain instances of'):
            Reader().run(filename)

    def test_validate(self):
        class TestModel(obj_model.Model):
            id = obj_model.StringAttribute(primary=True, unique=True)

        Model.Meta.attributes['test'] = obj_model.OneToOneAttribute(TestModel, related_name='a')
        with self.assertRaisesRegex(Exception, 'Relationships from `Model` not supported'):
            io.Writer.validate_implicit_relationships()
        Model.Meta.attributes.pop('test')

        Model.Meta.related_attributes['test'] = obj_model.OneToManyAttribute(TestModel, related_name='b')
        with self.assertRaisesRegex(Exception, 'Only one-to-one and many-to-one relationships are supported to `Model`'):
            io.Writer.validate_implicit_relationships()
        Model.Meta.related_attributes.pop('test')
