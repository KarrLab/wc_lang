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
    UnitDefinition, SBMLNamespaces, UNIT_KIND_SECOND, UNIT_KIND_MOLE)

from wc_lang.sbml.util import (wrap_libsbml, wrap_libsbml_pass_text, LibSBMLError, create_sbml_doc_w_fbc, wrap_libsbml_2,
    add_sbml_unit, create_sbml_parameter, init_sbml_model, SBML_LEVEL, SBML_VERSION, SBML_COMPATIBILITY_METHOD)


class TestSbml(unittest.TestCase):

    def setUp(self):
        # create an SBMLDocument that uses version 2 of the 'Flux Balance Constraints' extension
        self.document = create_sbml_doc_w_fbc()

    def test_SBML_wrap_libsbml(self):

        self.assertEqual(wrap_libsbml("LIBSBML_OPERATION_SUCCESS"), LIBSBML_OPERATION_SUCCESS)

        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml("1 +")
        self.assertIn("Syntax error in libsbml method", str(context.exception))

        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml("x")
        self.assertIn("NameError", str(context.exception))
        self.assertIn("'x' is not defined", str(context.exception))

        id = 'x'
        self.assertEqual(
            wrap_libsbml_pass_text("self.document.setIdAttribute", id), LIBSBML_OPERATION_SUCCESS)

        call = "self.document.setIdAttribute('..')"
        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml(call)
        self.assertIn('LibSBML returned error code', str(context.exception))
        self.assertIn("when executing '{}'".format(call), str(context.exception))

        call = "self.document.appendAnnotation(5)"
        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml(call)
        self.assertIn("in libsbml method call '{}'".format(call), str(context.exception))

        model = wrap_libsbml("self.document.createModel()")
        self.assertEqual(
            wrap_libsbml_pass_text("model.setTimeUnits", 'second'), LIBSBML_OPERATION_SUCCESS)

    def test_SBML_wrap_libsbml_2(self):

        id = 'x'
        self.assertEqual(
            wrap_libsbml_2(self.document.setIdAttribute, id), LIBSBML_OPERATION_SUCCESS)

        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml_2(self.document.setIdAttribute, '..')
        self.assertIn('LibSBML returned error code', str(context.exception))
        self.assertIn("when executing", str(context.exception))

        model = wrap_libsbml_2(self.document.createModel)
        self.assertEqual(
            wrap_libsbml_2(model.setTimeUnits, 'second'), LIBSBML_OPERATION_SUCCESS)

        self.assertEqual(
            wrap_libsbml_2(model.setTimeUnits, 'second',
                debug=True, returns_int=False, other=3), LIBSBML_OPERATION_SUCCESS)

    def test_init_sbml_model(self):
        sbml_model = init_sbml_model(self.document)

        # check the SBML document
        self.assertEqual(wrap_libsbml("self.document.checkConsistency()"), 0)
        self.assertEqual(wrap_libsbml("self.document.{}".format(SBML_COMPATIBILITY_METHOD)), 0)

        # check mmol_per_gDW_per_hr
        mmol_per_gDW_per_hr = wrap_libsbml("sbml_model.getUnitDefinition('mmol_per_gDW_per_hr')")
        printed_mmol_per_gDW_per_hr = wrap_libsbml("UnitDefinition.printUnits(mmol_per_gDW_per_hr)")
        compact_mmol_per_gDW_per_hr = wrap_libsbml("UnitDefinition.printUnits(mmol_per_gDW_per_hr, True)")
        self.assertIn('(0.001 mole)^1', compact_mmol_per_gDW_per_hr)
        self.assertIn('(3600 second)^-1', compact_mmol_per_gDW_per_hr)
        self.assertIn('(1 gram)^-1', compact_mmol_per_gDW_per_hr)

    def test_SBML_fbc(self):

        # create an SBMLDocument that uses version 2 of the 'Flux Balance Constraints' extension
        document = create_sbml_doc_w_fbc()

        id = 'x'
        self.assertEqual(
            wrap_libsbml_pass_text("document.setIdAttribute", id), LIBSBML_OPERATION_SUCCESS)


class TestLibsbmlInterface(unittest.TestCase):

    def setUp(self):
        sbmlns = wrap_libsbml('SBMLNamespaces({}, {}, "fbc", 2)'.format(SBML_LEVEL, SBML_VERSION))
        self.sbml_document = wrap_libsbml('SBMLDocument(sbmlns)')
        self.sbml_model = wrap_libsbml("self.sbml_document.createModel()")

        self.per_second_id = 'per_second'
        self.per_second = wrap_libsbml("self.sbml_model.createUnitDefinition()")
        wrap_libsbml_pass_text("self.per_second.setIdAttribute", self.per_second_id)
        add_sbml_unit(self.per_second, UNIT_KIND_SECOND, exponent=-1)

    def test_add_sbml_unit(self):
        per_second = wrap_libsbml("self.sbml_model.createUnitDefinition()")
        wrap_libsbml_pass_text("per_second.setIdAttribute", 'per_second')
        self.assertTrue(wrap_libsbml("per_second.hasRequiredAttributes()"))
        exp = -1
        default_scale=0
        default_multiplier=1.0
        unit = add_sbml_unit(per_second, UNIT_KIND_SECOND, exponent=exp)
        self.assertEqual(wrap_libsbml("unit.getExponent()", returns_int=True), exp)
        self.assertEqual(wrap_libsbml("unit.getKind()"), UNIT_KIND_SECOND)
        self.assertEqual(wrap_libsbml("unit.getScale()"), default_scale)
        self.assertEqual(wrap_libsbml("unit.getMultiplier()"), default_multiplier)

        strange_unit = wrap_libsbml("self.sbml_model.createUnitDefinition()")
        wrap_libsbml_pass_text("strange_unit.setIdAttribute", 'strange_unit')
        self.assertTrue(wrap_libsbml("strange_unit.hasRequiredAttributes()"))
        exp=-4; scale=3; mult=1.23
        unit = add_sbml_unit(strange_unit, UNIT_KIND_MOLE,
            exponent=exp, scale=scale, multiplier=mult)
        self.assertEqual(wrap_libsbml("unit.getExponent()", returns_int=True), exp)
        self.assertEqual(wrap_libsbml("unit.getKind()"), UNIT_KIND_MOLE)
        self.assertEqual(wrap_libsbml("unit.getScale()"), scale)
        self.assertEqual(wrap_libsbml("unit.getMultiplier()"), mult)

        with self.assertRaises(LibSBMLError) as context:
            unit = add_sbml_unit(strange_unit, -1)
        self.assertIn("LibSBML returned error code", str(context.exception))

    def test_create_sbml_parameter(self):
        id='id1'; name='name1'; value=13; constant=False
        parameter = create_sbml_parameter(self.sbml_model, id, name=name, value=value, constant=constant)
        self.assertTrue(wrap_libsbml("parameter.hasRequiredAttributes()"))
        self.assertEqual(wrap_libsbml("parameter.getIdAttribute()"), id)
        self.assertEqual(wrap_libsbml("parameter.getName()"), name)
        self.assertTrue(wrap_libsbml("parameter.isSetValue()"))
        self.assertFalse(wrap_libsbml("parameter.isSetUnits()"))
        self.assertEqual(wrap_libsbml("parameter.getValue()"), value)
        self.assertEqual(wrap_libsbml("parameter.getConstant()"), constant)

        # test defaults
        id = 'id2'
        parameter = create_sbml_parameter(self.sbml_model, id)
        self.assertEqual(wrap_libsbml("parameter.getIdAttribute()"), id)
        self.assertEqual(wrap_libsbml("parameter.getName()"), '')
        self.assertFalse(wrap_libsbml("parameter.isSetValue()"))
        self.assertEqual(wrap_libsbml("parameter.getConstant()"), True)

        # test units
        id = 'id3'
        parameter = create_sbml_parameter(self.sbml_model, id, units=self.per_second_id)
        self.assertTrue(wrap_libsbml("parameter.hasRequiredAttributes()"))
        self.assertTrue(wrap_libsbml("parameter.isSetUnits()"))
        self.assertEqual(wrap_libsbml("parameter.getUnits()"), self.per_second_id)

        # test Parameter id collision
        with self.assertRaises(ValueError) as context:
            parameter = create_sbml_parameter(self.sbml_model, id)
        self.assertIn("is already in use as a Parameter id", str(context.exception))

    def test_wrap_libsbml_pass_text(self):
        triple_double = '""" I can\'t mom """'
        triple_single = "''' hi mom '''"
        def f(a):
            return a
        self.assertEqual(wrap_libsbml_pass_text("f", triple_double), triple_double)
        self.assertEqual(wrap_libsbml_pass_text("f", triple_single), triple_single)

        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml_pass_text(None, 3)
        self.assertIn("3 isn't textual data", str(context.exception))
