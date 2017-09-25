""" Tests of SBML io

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import unittest
import os

from wc_lang.core import (Model, Taxon, Submodel, ObjectiveFunction, Compartment, SpeciesType,
    Species, Concentration, Reaction, ReactionParticipant, RateLaw, RateLawEquation,
    BiomassComponent, BiomassReaction, Parameter, Reference, CrossReference)

from wc_lang.sbml.io import Reader, Writer, SBMLExchange


class TestSbml(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_model.xlsx')

    '''
    def setUp(self):
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)
    '''
        
    def test_SBML_Exchange(self):
        document = SBMLExchange.write([], [Compartment])

