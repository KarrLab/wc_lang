""" Tests of SBML

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import unittest
import os

# "from libsbml import *" generates "NameError: Unknown C global variable" in pytest,
# presumably from the SWIG wrapper: http://web.mit.edu/svn/src/swig-1.3.25/Lib/python/pyinit.swg
from libsbml import (LIBSBML_OPERATION_SUCCESS, SBMLDocument, OperationReturnValue_toString)

from wc_lang.core import (Model, Taxon, Submodel, ObjectiveFunction, Compartment, SpeciesType,
    Species, Concentration, Reaction, ReactionParticipant, RateLaw, RateLawEquation,
    BiomassComponent, BiomassReaction, Parameter, Reference, CrossReference)

from wc_lang.sbml import wrap_libsbml, LibSBMLError, Reader, Writer, SBMLExchange


class TestSbml(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_model.xlsx')

    '''
    def setUp(self):
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)
    '''
        
    def test_SBML_wrap_libsbml(self):

        self.assertEqual(wrap_libsbml("LIBSBML_OPERATION_SUCCESS"), LIBSBML_OPERATION_SUCCESS)

        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml("1 +")
        self.assertIn("Syntax error in libsbml method", str(context.exception))

        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml("x")
        self.assertIn("NameError", str(context.exception))
        self.assertIn("'x' is not defined", str(context.exception))

        try:
            document = SBMLDocument(3, 1)
        except ValueError:
            raise SystemExit("'SBMLDocument(3, 1)' fails")
        
        id = 'x'
        self.assertEqual(
            wrap_libsbml("document.setIdAttribute('{}')".format(id)), LIBSBML_OPERATION_SUCCESS)

        id = '..'
        call = "document.setIdAttribute('{}')".format(id)
        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml(call)
        self.assertIn('LibSBML returned error code', str(context.exception))
        self.assertIn("when executing '{}'".format(call), str(context.exception))

        call = "document.appendAnnotation(5)"
        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml(call)
        self.assertIn("in libsbml method call '{}'".format(call), str(context.exception))

        model = wrap_libsbml("document.createModel()")
        self.assertEqual(wrap_libsbml("model.setTimeUnits('second')"), LIBSBML_OPERATION_SUCCESS)

    def test_SBML_Exchange(self):
        # document = SBMLExchange.write([self.model.get_compartments()], [Compartment])
        document = SBMLExchange.write([], [Compartment])
