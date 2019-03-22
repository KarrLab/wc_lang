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
from wc_lang.core import Model, WcLangWarning
from wc_lang.io import Reader
from wc_lang.sbml import io as sbml_io
from wc_lang.sbml import util as sbml_util
from wc_lang.transform.prep_for_sbml import PrepForSbmlTransform
from wc_utils.util.ontology import wcm_ontology
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
        submodel = model.submodels.create(id='Metabolism', framework=wcm_ontology['WCM:dynamic_flux_balance_analysis'])
        model.reactions.create(id='rxn_1', submodel=submodel, flux_bound_units=unit_registry.parse_units('M s^-1'))
        model.reactions.create(id='rxn_1', submodel=submodel, flux_bound_units=unit_registry.parse_units('M s^-1'))

        with self.assertRaisesRegex(sbml_util.LibSbmlError, 'Document is invalid'):
            with self.assertRaisesRegex(WcLangWarning, 'Model is invalid'):
                sbml_io.SbmlExporter.run(model)

    def test_write_read(self):
        submodels = sbml_io.SbmlWriter().run(self.model, self.dirname)
        model = sbml_io.SbmlReader().run(self.dirname)

        sbml_compat_model = PrepForSbmlTransform().run(self.model.copy())
        print(model.difference(sbml_compat_model))
        self.assertTrue(model.is_equal(sbml_compat_model))

    def test_writer_errors(self):
        with mock.patch('libsbml.writeSBMLToFile', return_value=False):
            with self.assertRaisesRegex(ValueError, ' could not be written to '):
                sbml_io.SbmlWriter().run(self.model, self.dirname)
