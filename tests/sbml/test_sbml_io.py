""" Tests of SBML io

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

from math import isnan
from six import iteritems
import mock
import os
import shutil
import tempfile
import unittest

import libsbml
from libsbml import Compartment as libsbmlCompartment
from libsbml import Species as libsbmlSpecies
from libsbml import Parameter as libsbmlParameter
from libsbml import Reaction as libsbmlReaction
from libsbml import Model as libsbmlModel
from libsbml import Objective as libsbmlObjective

from obj_model.utils import get_component_by_id
from wc_lang.core import (SubmodelAlgorithm, Model, Taxon, Submodel, ObjectiveFunction, Compartment,
                          Species, Concentration, Reaction, ReactionParticipant, RateLaw, RateLawEquation,
                          BiomassComponent, BiomassReaction, Parameter, Reference, CrossReference)
from wc_lang.prepare import PrepareModel, CheckModel

from wc_lang.sbml.util import wrap_libsbml, get_SBML_compatibility_method
from wc_lang.io import Reader
import wc_lang.sbml.io as sbml_io

# ignore 'setting concentration' warnings
import warnings
warnings.filterwarnings('ignore', '.*setting concentration.*', )


def check_document_against_model(sbml_document, wc_lang_model, test_case):
    """ Compare an SBML document against a wc lang model.

    Check that the `sbml_document` is consistent with the `wc_lang_model`. Many classes and some
    instances in `wc_lang_model` may be missing from `sbml_document`.

    Args:
        sbml_document (:obj:`libsbml.SBMLDocument`): a libsbml SBMLDocument
        wc_lang_model (:obj:`Model`): a wc lang `Model` with species, reactions, compartments or parameters
        test_case (:obj:`unittest.TestCase`): a unittest TestCase
    """
    for element in sbml_document.getListOfAllElements():

        # compartments
        if isinstance(element, libsbmlCompartment):
            wc_lang_compartment = get_component_by_id(wc_lang_model.get_compartments(),
                                                      element.getIdAttribute())
            test_case.assertEqual(element.getName(), wc_lang_compartment.name)
            test_case.assertEqual(element.getSpatialDimensions(), 3)
            test_case.assertEqual(element.getSize(), wc_lang_compartment.initial_volume)
            continue
            # not checking: comments

        # parameters
        if isinstance(element, libsbmlParameter):
            prefix = 'parameter_'
            if element.getIdAttribute().startswith(prefix):
                # only parameters that start with 'parameter' are wc_lang Parameters
                wc_lang_id = element.getIdAttribute()[len(prefix):]
                wc_lang_parameter = get_component_by_id(wc_lang_model.get_parameters(), wc_lang_id)
                test_case.assertEqual(element.getName(), wc_lang_parameter.name)
                test_case.assertEqual(element.getValue(), wc_lang_parameter.value)
            continue
            # not checking: units, SBML parameters not in wc_lang Parameters

        if isinstance(element, libsbmlSpecies):
            # because Species.id() is a method, get_component_by_id() cannot search it
            species_index = {species.id(): species for species in wc_lang_model.get_species()}
            wc_lang_id = Species.xml_id_to_id(element.getIdAttribute())
            if not wc_lang_id in species_index:
                test_case.assertFail("Cannot find species id '{}' in species index".format(wc_lang_id))
            wc_lang_species = species_index[wc_lang_id]
            test_case.assertEqual(element.getName(), wc_lang_species.species_type.name)
            test_case.assertEqual(element.getCompartment(), wc_lang_species.compartment.id)
            test_case.assertEqual(element.getInitialConcentration(), wc_lang_species.concentration.value)
            continue
            # not checking: comments

        if isinstance(element, libsbmlReaction):
            wc_lang_reaction = get_component_by_id(wc_lang_model.get_reactions(), element.getIdAttribute())
            # test Reaction
            if wc_lang_reaction:
                test_case.assertEqual(element.getName(), wc_lang_reaction.name)
                test_case.assertEqual(element.getReversible(), wc_lang_reaction.reversible)
                test_case.assertEqual(element.getFast(), False)
                test_case.assertEqual(element.getCompartment(), wc_lang_reaction.submodel.compartment.id)
                continue
                # not checking: participants and flux bounds

            wc_lang_biomass_reaction = get_component_by_id(wc_lang_model.get_biomass_reactions(),
                                                           element.getIdAttribute())
            # test BiomassReaction
            if wc_lang_biomass_reaction:
                test_case.assertEqual(element.getName(), wc_lang_biomass_reaction.name)
                test_case.assertEqual(element.getReversible(), False)
                test_case.assertEqual(element.getFast(), False)
                test_case.assertEqual(element.getCompartment(), wc_lang_biomass_reaction.compartment.id)
                continue
                # not checking: components, flux bounds, and comments

        if isinstance(element, libsbmlModel):
            # test a submodel
            wc_lang_submodel = get_component_by_id(wc_lang_model.get_submodels(), element.getIdAttribute())
            if wc_lang_submodel.name:
                test_case.assertEqual(element.getName(), wc_lang_submodel.name)
            continue
            # not checking: comments

        if isinstance(element, libsbmlObjective):
            # test an ObjectiveFunction
            test_case.assertEqual(element.getType(), 'maximize')
            test_case.assertEqual(element.getIdAttribute(), ObjectiveFunction.ACTIVE_OBJECTIVE)
            continue
            # not checking: reactions, or biomass_reactions

            # TODO: check remaining elements


class TestSbml(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), '..', 'fixtures', 'example-model.xlsx')

    def setUp(self):
        self.dirname = tempfile.mkdtemp()
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)
        PrepareModel(self.model).run()
        CheckModel(self.model).run()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def check_sbml_doc(self, sbml_doc):

        # if checkConsistency() returns some errors, print them
        for i in range(wrap_libsbml(sbml_doc.checkConsistency, returns_int=True)):
            print(sbml_doc.getError(i).getShortMessage())
            print(sbml_doc.getError(i).getMessage())
        self.assertEqual(wrap_libsbml(sbml_doc.checkConsistency, returns_int=True), 0)

        # if 0<getNumErrors, print them
        for i in range(wrap_libsbml(sbml_doc.getNumErrors, returns_int=True)):
            print(sbml_doc.getError(i).getShortMessage())
            print(sbml_doc.getError(i).getMessage())
        self.assertEqual(wrap_libsbml(sbml_doc.getNumErrors, returns_int=True), 0)

    def test_SBML_Exchange(self):
        for submodel in self.model.get_submodels():
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                sbml_document = sbml_io.SBMLExchange.write_submodel(submodel)

                self.assertEqual(wrap_libsbml(get_SBML_compatibility_method(sbml_document),
                                              returns_int=True), 0)
                self.check_sbml_doc(sbml_document)
                check_document_against_model(sbml_document, self.model, self)

    def test_SBML_Exchange_warning(self):
        model = Model(id='model')
        met_submodel = model.submodels.create(id='Metabolism', algorithm=SubmodelAlgorithm.ssa)
        met_submodel.parameters.append(Parameter(id='param_1'))
        met_submodel.parameters.append(Parameter(id='param_1'))

        with self.assertRaisesRegexp(UserWarning, 'Some data will not be written because objects are not valid'):
            warnings.simplefilter("ignore")
            warnings.simplefilter("error", UserWarning)
            sbml_document = sbml_io.SBMLExchange.write_submodel(met_submodel)
            warnings.resetwarnings()

    def test_writer(self):
        for algorithms in [None, [SubmodelAlgorithm.dfba]]:
            sbml_documents = sbml_io.Writer.run(self.model, algorithms=algorithms)
            try:
                paths = sbml_io.Writer.run(self.model, algorithms=algorithms, path=self.dirname)
            except Exception as e:
                self.fail("Unexpected sbml_io.Writer.run() exception '{}'".format(e))
            for submodel_id, path in zip(sbml_documents.keys(), paths):

                document = libsbml.SBMLReader().readSBML(path)
                self.check_sbml_doc(document)

                self.assertEqual(document.toSBML(), sbml_documents[submodel_id].toSBML())
                check_document_against_model(document, self.model, self)

    def test_writer_errors(self):
        root_path = os.path.join(os.path.dirname(__file__), 'no_such_dir', 'example-model')
        with self.assertRaises(ValueError) as context:
            sbml_io.Writer.run(self.model, path=root_path)
        self.assertIn('cannot write to directory', str(context.exception))
        self.assertIn('no_such_dir', str(context.exception))

        with self.assertRaisesRegexp(ValueError, 'No submodel.algorithm in algorithms'):
            sbml_io.Writer.run(self.model, algorithms=[])

        with mock.patch('libsbml.writeSBMLToFile', return_value=False):
            with self.assertRaisesRegexp(ValueError, ' could not be written to '):
                sbml_io.Writer.run(self.model, path=self.dirname)
