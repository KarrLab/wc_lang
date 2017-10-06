""" Tests of SBML io

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import unittest
import os
from math import isnan
from libsbml import readSBMLFromString, writeSBMLToFile

from wc_lang.core import (SubmodelAlgorithm, Model, Taxon, Submodel, ObjectiveFunction, Compartment,
    Species, Concentration, Reaction, ReactionParticipant, RateLaw, RateLawEquation,
    BiomassComponent, BiomassReaction, Parameter, Reference, CrossReference)
from wc_lang.sbml.util import wrap_libsbml, SBML_COMPATIBILITY_METHOD
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
                print("setting concentration for {} to 0.0".format(specie.id()))
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
                    objects + [submodel, submodel.get_objective_function()])

        # TODO: either compile & use the libsbml source trunk, or install the next release of libsbml
        self.assertEqual(wrap_libsbml("sbml_document.{}".format(SBML_COMPATIBILITY_METHOD)), 0)
        workaround_document = wrap_libsbml("readSBMLFromString(sbml_document.toSBML())")
        self.assertEqual(wrap_libsbml("workaround_document.checkConsistency()"), 0)
        for i in range(wrap_libsbml("workaround_document.checkConsistency()")):
            print(workaround_document.getError(i).getShortMessage())
            print(workaround_document.getError(i).getMessage())
