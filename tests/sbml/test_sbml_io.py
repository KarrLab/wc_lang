""" Tests of SBML io

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import unittest
import os
from math import isnan
from six import iteritems

from libsbml import readSBMLFromString, writeSBMLToFile, SBMLReader, SBMLDocument
from libsbml import Compartment as libsbmlCompartment
from libsbml import Species as libsbmlSpecies
from libsbml import Parameter as libsbmlParameter
from libsbml import Reaction as libsbmlReaction

from obj_model.utils import get_component_by_id
from wc_lang.core import (SubmodelAlgorithm, Model, Taxon, Submodel, ObjectiveFunction, Compartment,
    Species, Concentration, Reaction, ReactionParticipant, RateLaw, RateLawEquation,
    BiomassComponent, BiomassReaction, Parameter, Reference, CrossReference)
from wc_lang.sbml.util import wrap_libsbml, get_SBML_compatibility_method
from wc_lang.io import Reader
import wc_lang.sbml.io as sbml_io

def check_document_against_model(sbml_document, wc_lang_model, test_case):
    """ Compare an SBML document against a wc lang model.

    Check that the `sbml_document` is consistent with the `wc_lang_model`. Many classes and some
    instances in `wc_lang_model` may be missing from `sbml_document`.

    Args:
        sbml_document (:obj:`SBMLDocument`): a libsbml SBMLDocument
        wc_lang_model (:obj:`Model`): a wc lang `Model` with species, reactions, compartments or parameters
        test_case (:obj:`unittest.TestCase`): a TestCase
    """
    for element in sbml_document.getListOfAllElements():
        # compartments
        if isinstance(element, libsbmlCompartment):
            wc_lang_compartment = get_component_by_id(wc_lang_model.get_compartments(),
                element.getIdAttribute())
            test_case.assertEqual(element.getName(), wc_lang_compartment.name)
            test_case.assertEqual(element.getSpatialDimensions(), 3)
            test_case.assertEqual(element.getSize(), wc_lang_compartment.initial_volume)
            # not checking: comments

        # parameters
        elif isinstance(element, libsbmlParameter):
            prefix = 'parameter_'
            if element.getIdAttribute().startswith(prefix):
                # only parameters that start with 'parameter' are wc_lang Parameters
                wc_lang_id = element.getIdAttribute()[len(prefix):]
                wc_lang_parameter = get_component_by_id(wc_lang_model.get_parameters(), wc_lang_id)
                test_case.assertEqual(element.getName(), wc_lang_parameter.name)
                test_case.assertEqual(element.getValue(), wc_lang_parameter.value)
            # not checking: units, SBML parameters not in wc_lang Parameters

        elif isinstance(element, libsbmlSpecies):
            # because Species.id() is a method, get_component_by_id() cannot search it
            species_index = {species.id(): species for species in wc_lang_model.get_species()}
            wc_lang_id = Species.xml_id_to_id(element.getIdAttribute())
            if not wc_lang_id in species_index:
                test_case.assertFail("Cannot find species id '{}' in species index".format(wc_lang_id))
            wc_lang_species = species_index[wc_lang_id]
            test_case.assertEqual(element.getName(), wc_lang_species.species_type.name)
            test_case.assertEqual(element.getCompartment(), wc_lang_species.compartment.id)
            test_case.assertEqual(element.getInitialConcentration(), wc_lang_species.concentration.value)
            # not checking: comments

        elif isinstance(element, libsbmlReaction):
            wc_lang_reaction = get_component_by_id(wc_lang_model.get_reactions(), element.getIdAttribute())
            test_case.assertEqual(element.getName(), wc_lang_reaction.name)
            test_case.assertEqual(element.getReversible(), wc_lang_reaction.reversible)
            test_case.assertEqual(element.getFast(), False)
            test_case.assertEqual(element.getCompartment(), wc_lang_reaction.submodel.compartment.id)
            # not checking: participants and flux bounds

            # TODO: check remaining elements, and ObjectiveFunction

class TestSbml(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), '..', 'fixtures', 'example-model.xlsx')

    def setUp(self):
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)

        # hack in concentration of 0 until we have real consistency checking
        # TODO: replace with real consistency checking
        for specie in self.model.get_species():
            if specie.concentration is None:
                # TODO: make this a warning
                # print("setting concentration for {} to 0.0".format(specie.id()))
                specie.concentrations = Concentration(species=specie, value=0.0)

        # TODO: replace with real code for setting bounds
        default_min_flux_bound = 0
        default_max_flux_bound = 1000
        for rxn in self.model.get_reactions():
            if isnan(rxn.min_flux):
                if rxn.reversible:
                    rxn.min_flux = -default_max_flux_bound
                else:
                    rxn.min_flux = default_min_flux_bound
            if isnan(rxn.max_flux):
                rxn.max_flux = default_max_flux_bound

    def test_SBML_Exchange(self):
        objects = \
            self.model.get_compartments() + \
            self.model.get_species() + \
            self.model.get_reactions() + \
            self.model.get_parameters()
        for submodel in self.model.get_submodels():
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                sbml_document = sbml_io.SBMLExchange.write(
                    objects + [submodel, submodel.objective_function])

        # TODO: avoid workaround by installing libsbml>15.5.0
        self.assertEqual(wrap_libsbml(get_SBML_compatibility_method(sbml_document)), 0)
        sbml_string = wrap_libsbml(sbml_document.toSBML)
        workaround_document = wrap_libsbml(readSBMLFromString, sbml_string)
        # if checkConsistency() returns some errors, print them
        for i in range(wrap_libsbml(workaround_document.checkConsistency)):
            print(workaround_document.getError(i).getShortMessage())
            print(workaround_document.getError(i).getMessage())
        self.assertEqual(wrap_libsbml(workaround_document.checkConsistency), 0)
        check_document_against_model(sbml_document, self.model, self)

    def test_writer(self):
        root_path = os.path.join(os.path.dirname(__file__), 'fixtures', 'example-model')
        for algorithms in [None, [SubmodelAlgorithm.dfba]]:
            sbml_documents = sbml_io.Writer.run(self.model, algorithms=algorithms, path=None)
            try:
                paths = sbml_io.Writer.run(self.model, algorithms=algorithms, path=root_path)
            except Exception as e:
                self.fail("Unexpected sbml_io.Writer.run() exception '{}'".format(e))
            for submodel_id,path in zip(sbml_documents.keys(), paths):
                document = SBMLReader().readSBML(path)
                for i in range(wrap_libsbml(document.checkConsistency)):
                    print(document.getError(i).getShortMessage())
                    print(document.getError(i).getMessage())
                self.assertEqual(wrap_libsbml(document.checkConsistency), 0)
                for i in range(wrap_libsbml(document.getNumErrors)):
                    print(document.getError(i).getShortMessage())
                    print(document.getError(i).getMessage())
                self.assertEqual(wrap_libsbml(document.getNumErrors), 0)
                self.assertEqual(document.toSBML(), sbml_documents[submodel_id].toSBML())
                check_document_against_model(document, self.model, self)

    def test_writer_errors(self):
        root_path = os.path.join(os.path.dirname(__file__), 'no_such_dir', 'example-model')
        with self.assertRaises(ValueError) as context:
            sbml_io.Writer.run(self.model, path=root_path)
        self.assertIn('cannot write to directory', str(context.exception))
        self.assertIn('no_such_dir', str(context.exception))
