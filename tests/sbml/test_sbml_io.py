""" Tests of SBML io

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2019-03-21
:Copyright: 2017-2019, Karr Lab
:License: MIT
"""

import libsbml
import mock
import os
import shutil
import tempfile
import unittest
from wc_lang.core import (Model, Submodel, Compartment, Reaction, FluxBounds, RateLaw, RateLawDirection,
                          DfbaObjective, DfbaObjectiveExpression, WcLangWarning, InitVolume,
                          ChemicalStructure, ChemicalStructureFormat, ChemicalStructureAlphabet)
from wc_lang.io import Reader
from wc_lang.sbml import io as sbml_io
from wc_lang.sbml import util as sbml_util
from wc_lang.transform.prep_for_sbml import PrepForSbmlTransform
from wc_onto import onto
from wc_utils.util.units import unit_registry


class SbmlIoTestCase(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), '..', 'fixtures', 'sbml-io.xlsx')

    def setUp(self):
        self.dirname = tempfile.mkdtemp()
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)[Model][0]

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_export_import(self):
        model = PrepForSbmlTransform().run(self.model)
        core, submodels = model.submodels.gen_models()

        all_submodels = [core] + submodels
        all_submodels_2 = []
        for i_submodel, submodel in enumerate(all_submodels):

            sbml_doc = sbml_io.SbmlExporter.run(submodel)
            submodel_2 = sbml_io.SbmlImporter.run(sbml_doc)
            self.assertTrue(submodel_2.is_equal(submodel))
            all_submodels_2.append(submodel_2)

        model_2 = all_submodels_2[0]
        for submodel_2 in all_submodels_2[1:]:
            model_2.merge(submodel_2)
        self.assertTrue(model_2.is_equal(model))

    def test_SbmlExporter_error(self):
        model = Model(id='model')
        model.submodels.create(id='Metabolism', framework=onto['WC:dynamic_flux_balance_analysis'])
        model.submodels.create(id='Metabolism_2', framework=onto['WC:dynamic_flux_balance_analysis'])
        with self.assertRaisesRegex(ValueError, 'Only 1 submodel can be encoded'):
            sbml_io.SbmlExporter.run(model)

        model = Model(id='model')
        submodel = model.submodels.create(id='Metabolism', framework=onto['WC:dynamic_flux_balance_analysis'])
        model.reactions.create(id='rxn_1', submodel=submodel, flux_bounds=FluxBounds(units=unit_registry.parse_units('M s^-1')))
        model.reactions.create(id='rxn_1', submodel=submodel, flux_bounds=FluxBounds(units=unit_registry.parse_units('M s^-1')))
        with self.assertRaisesRegex(sbml_util.LibSbmlError, 'Document is invalid'):
            with self.assertWarnsRegex(WcLangWarning, 'Model is invalid'):
                sbml_io.SbmlExporter.run(model)

    def test_write_read(self):
        shutil.rmtree(self.dirname)

        sbml_io.SbmlWriter().run(self.model, self.dirname)
        model = sbml_io.SbmlReader().run(self.dirname)

        sbml_compat_model = PrepForSbmlTransform().run(self.model.copy())
        self.assertTrue(model.is_equal(sbml_compat_model))

    def test_writer_errors(self):
        model = Model(version=None)
        with self.assertWarnsRegex(WcLangWarning, 'Model is invalid:'):
            sbml_io.SbmlWriter().run(model, self.dirname)

        with mock.patch('libsbml.writeSBMLToFile', return_value=False):
            with self.assertRaisesRegex(ValueError, ' could not be written to '):
                sbml_io.SbmlWriter().run(self.model, self.dirname)


class SbmlIoInCoreTestCase(unittest.TestCase):
    def test_Compartment_to_from_sbml(self):
        c = Compartment(id='c', name='cytosol', geometry=onto['WC:3D_compartment'])

        sbml_doc = sbml_util.LibSbmlInterface.create_doc()
        sbml_model = sbml_util.LibSbmlInterface.create_model(sbml_doc)

        sbml_c = c.export_to_sbml(sbml_model)

        c_2 = Compartment()
        c_2.import_from_sbml(sbml_c)

        self.assertTrue(c_2.is_equal(c))

    def test_Compartment_to_sbml_error(self):
        c = Compartment(id='c', name='cytosol', geometry=None)
        sbml_doc = sbml_util.LibSbmlInterface.create_doc()
        sbml_model = sbml_util.LibSbmlInterface.create_model(sbml_doc)
        with self.assertRaisesRegex(ValueError, 'Unsupported geometry'):
            c.export_to_sbml(sbml_model)

        c = Compartment(id='c', name='cytosol', geometry=onto['WC:3D_compartment'])
        sbml_doc = sbml_util.LibSbmlInterface.create_doc()
        sbml_model = sbml_util.LibSbmlInterface.create_model(sbml_doc)
        sbml_c = c.export_to_sbml(sbml_model)
        sbml_util.LibSbmlInterface.call_libsbml(sbml_c.setSpatialDimensions, 2)
        c_2 = Compartment()
        with self.assertRaisesRegex(ValueError, 'Unsupported spatial dimensions'):
            c_2.import_from_sbml(sbml_c)

    def test_Compartment_from_sbml(self):
        model = Model()
        init_volume = InitVolume(mean=1.2)
        c_1 = model.compartments.create(id='c_1', name='cytosol', geometry=onto['WC:3D_compartment'], init_volume=init_volume)
        c_2 = model.compartments.create(id='c_2', name='cytosol', geometry=onto['WC:3D_compartment'], init_volume=init_volume)

        sbml_doc = sbml_util.LibSbmlInterface.create_doc()
        sbml_model = sbml_util.LibSbmlInterface.create_model(sbml_doc)

        sbml_c_1 = c_1.export_to_sbml(sbml_model)
        sbml_c_2 = c_2.export_to_sbml(sbml_model)

        model_2 = Model()
        c_1_b = model_2.compartments.create()
        c_1_b.import_from_sbml(sbml_c_1)
        c_2_b = model_2.compartments.create()
        c_2_b.import_from_sbml(sbml_c_2)

        self.assertTrue(c_1_b.is_equal(c_1))
        self.assertTrue(c_2_b.is_equal(c_2))

    def test_Species_from_sbml(self):
        model = Model(id='mdl', created=None, updated=None)

        c = model.compartments.create(id='c')

        st_1 = model.species_types.create(id='st_1', structure=ChemicalStructure(
            value='AAA', format=ChemicalStructureFormat.BpForms, alphabet=ChemicalStructureAlphabet.dna))
        st_1.structure.empirical_formula = st_1.structure.get_structure().get_formula()
        st_1.structure.molecular_weight = st_1.structure.get_structure().get_mol_wt()
        st_1.structure.charge = st_1.structure.get_structure().get_charge()
        st_2 = model.species_types.create(id='st_2')

        s_1 = model.species.create(species_type=st_1, compartment=c)
        s_2 = model.species.create(species_type=st_2, compartment=c)
        s_1.id = s_1.gen_id()
        s_2.id = s_2.gen_id()

        sbml_doc = sbml_util.LibSbmlInterface.create_doc()
        sbml_model = sbml_util.LibSbmlInterface.create_model(sbml_doc)

        s_1_sbml = s_1.export_to_sbml(sbml_model)
        s_2_sbml = s_2.export_to_sbml(sbml_model)
        s_1.export_relations_to_sbml(sbml_model, s_1_sbml)
        s_2.export_relations_to_sbml(sbml_model, s_2_sbml)

        model_2 = Model(id='mdl', created=None, updated=None)
        c_b = model_2.compartments.create(id='c')
        s_1_b = model_2.species.create()
        s_2_b = model_2.species.create()
        s_1_b.import_from_sbml(s_1_sbml)
        s_2_b.import_from_sbml(s_2_sbml)
        objs = {Compartment: {c_b.id: c_b}}
        s_1_b.import_relations_from_sbml(s_1_sbml, objs=objs)
        s_2_b.import_relations_from_sbml(s_2_sbml, objs=objs)

        self.assertTrue(s_1_b.is_equal(s_1))
        self.assertTrue(s_2_b.is_equal(s_2))
        self.assertTrue(model_2.is_equal(model))

    def test_Reaction_to_sbml(self):
        rxn = Reaction(submodel=Submodel(), reversible=True, rate_laws=[RateLaw(direction=RateLawDirection.backward)])

        sbml_doc = sbml_util.LibSbmlInterface.create_doc()
        sbml_model = sbml_util.LibSbmlInterface.create_model(sbml_doc)

        with self.assertRaisesRegex(ValueError, 'must be split'):
            rxn.export_to_sbml(sbml_model)

    def test_RateLaw_from_sbml(self):
        self.assertEqual(RateLaw._import_relations_from_sbml_units_transform('k_cat'), 'k_cat')
        self.assertEqual(RateLaw._import_relations_from_sbml_units_transform('k_cat * 1 mole'), 'k_cat')
        self.assertEqual(RateLaw._import_relations_from_sbml_units_transform('(k_cat) * 1 mole'), 'k_cat')

    def test_DfbaObjective_to_sbml_warning(self):
        obj = DfbaObjective(submodel=Submodel(), expression=DfbaObjectiveExpression())
        obj.expression._parsed_expression = mock.Mock(is_linear=False)

        with self.assertWarnsRegex(WcLangWarning, "doesn't support the non-linear objective"):
            obj.export_to_sbml(None)
