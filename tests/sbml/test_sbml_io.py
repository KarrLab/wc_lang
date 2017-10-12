""" Tests of SBML io

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import unittest
import os
from math import isnan
from libsbml import readSBMLFromString, writeSBMLToFile, SBMLReader
from six import iteritems

from wc_lang.core import (SubmodelAlgorithm, Model, Taxon, Submodel, ObjectiveFunction, Compartment,
    Species, Concentration, Reaction, ReactionParticipant, RateLaw, RateLawEquation,
    BiomassComponent, BiomassReaction, Parameter, Reference, CrossReference)
from wc_lang.sbml.util import wrap_libsbml, get_SBML_compatibility_method
from wc_lang.io import Reader
import wc_lang.sbml.io as sbml_io

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
        for i in range(wrap_libsbml(workaround_document.checkConsistency)):
            print(workaround_document.getError(i).getShortMessage())
            print(workaround_document.getError(i).getMessage())
        self.assertEqual(wrap_libsbml(workaround_document.checkConsistency), 0)

        # check some data in sbml document
        # TODO: check one instance of each class

    def test_writer(self):
        root_path = os.path.join(os.path.dirname(__file__), 'fixtures', 'example-model')
        # for algorithms in [None, [SubmodelAlgorithm.dfba, SubmodelAlgorithm.ssa]]:
        for algorithms in [None]:
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


    def test_writer_errors(self):
        root_path = os.path.join(os.path.dirname(__file__), 'no_such_dir', 'example-model')
        with self.assertRaises(ValueError) as context:
            sbml_io.Writer.run(self.model, path=root_path)
        self.assertIn('cannot write to directory', str(context.exception))
        self.assertIn('no_such_dir', str(context.exception))
