""" Tests of SBML utils

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import unittest
import os

# "from libsbml import *" generates "NameError: Unknown C global variable" in pytest,
# presumably from the SWIG wrapper: http://web.mit.edu/svn/src/swig-1.3.25/Lib/python/pyinit.swg
from libsbml import (LIBSBML_OPERATION_SUCCESS, SBMLDocument, OperationReturnValue_toString,
    SBMLNamespaces)

from wc_lang.sbml.util import wrap_libsbml, LibSBMLError


class TestSbml(unittest.TestCase):

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

        call = "document.setIdAttribute('..')"
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

    def test_init_sbml_model(self):
        # TODO
        pass

    def test_SBML_fbc(self):

        try:
            # use uses the SBML Level 3 Flux Balance Constraints package
            sbmlns = SBMLNamespaces(3, 2, "fbc", 2);
            document = SBMLDocument(sbmlns);
            # mark the fbc package required
            document.setPackageRequired("fbc", True)
        except ValueError:
            raise SystemExit("'SBMLNamespaces(3, 2, 'fbc', 2) fails")

        id = 'x'
        self.assertEqual(
            wrap_libsbml("document.setIdAttribute('{}')".format(id)), LIBSBML_OPERATION_SUCCESS)
