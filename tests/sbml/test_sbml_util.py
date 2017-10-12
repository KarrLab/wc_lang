""" Tests of SBML utils

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import unittest
import os
import six

# "from libsbml import *" generates "NameError: Unknown C global variable" in pytest,
# presumably from the SWIG wrapper: http://web.mit.edu/svn/src/swig-1.3.25/Lib/python/pyinit.swg
from libsbml import (LIBSBML_OPERATION_SUCCESS, SBMLDocument, OperationReturnValue_toString,
    UnitDefinition, SBMLNamespaces, UNIT_KIND_SECOND, UNIT_KIND_MOLE)

from wc_lang.sbml.util import (wrap_libsbml, LibSBMLError, create_sbml_doc_w_fbc, add_sbml_unit,
    create_sbml_parameter, init_sbml_model, SBML_LEVEL, SBML_VERSION, get_SBML_compatibility_method)


class TestSbml(unittest.TestCase):

    def setUp(self):
        # create an SBMLDocument that uses version 2 of the 'Flux Balance Constraints' extension
        self.document = create_sbml_doc_w_fbc()

    def test_SBML_wrap_libsbml(self):

        # id = 'test_id'
        id = u'test_id'
        self.assertEqual(
            wrap_libsbml(self.document.setIdAttribute, id), LIBSBML_OPERATION_SUCCESS)
        self.assertEqual(
            str(wrap_libsbml(self.document.getIdAttribute)), str(id))

        model = wrap_libsbml(self.document.createModel)
        self.assertEqual(
            wrap_libsbml(model.setTimeUnits, 'second'), LIBSBML_OPERATION_SUCCESS)

        self.assertEqual(
            wrap_libsbml(model.setTimeUnits, 'second',
                debug=True, returns_int=False, other=3), LIBSBML_OPERATION_SUCCESS)

        self.assertEqual(
            wrap_libsbml(model.setTimeUnits, 'second',
                debug=True, returns_int=False, other=3), LIBSBML_OPERATION_SUCCESS)

        self.assertEqual(
            wrap_libsbml(self.document.getNumErrors, returns_int=False), 0)

        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml(self.document.getNumErrors, 'no arg')
        self.assertIn('Error', str(context.exception))
        self.assertIn("in libsbml method call", str(context.exception))

        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml(self.document.setIdAttribute, '..')
        self.assertIn('LibSBML returned error code', str(context.exception))
        self.assertIn("when executing", str(context.exception))

        with self.assertRaises(LibSBMLError) as context:
            wrap_libsbml(self.document.getAnnotation)
        self.assertIn('libsbml returned None when executing', str(context.exception))

    def test_init_sbml_model(self):
        sbml_model = init_sbml_model(self.document)

        # check the SBML document
        self.assertEqual(wrap_libsbml(self.document.checkConsistency), 0)
        self.assertEqual(wrap_libsbml(get_SBML_compatibility_method(self.document)), 0)

        # check mmol_per_gDW_per_hr
        mmol_per_gDW_per_hr = wrap_libsbml(sbml_model.getUnitDefinition, 'mmol_per_gDW_per_hr')
        printed_mmol_per_gDW_per_hr = wrap_libsbml(UnitDefinition.printUnits, mmol_per_gDW_per_hr)
        compact_mmol_per_gDW_per_hr = wrap_libsbml(UnitDefinition.printUnits, mmol_per_gDW_per_hr, True)
        self.assertIn('(0.001 mole)^1', compact_mmol_per_gDW_per_hr)
        self.assertIn('(3600 second)^-1', compact_mmol_per_gDW_per_hr)
        self.assertIn('(1 gram)^-1', compact_mmol_per_gDW_per_hr)

    def test_SBML_fbc(self):

        # create an SBMLDocument that uses version 2 of the 'Flux Balance Constraints' extension
        document = create_sbml_doc_w_fbc()

        id = 'x'
        self.assertEqual(
            wrap_libsbml(document.setIdAttribute, id), LIBSBML_OPERATION_SUCCESS)


class TestLibsbmlInterface(unittest.TestCase):

    def setUp(self):
        sbmlns = wrap_libsbml(SBMLNamespaces, SBML_LEVEL, SBML_VERSION, "fbc", 2)
        self.sbml_document = wrap_libsbml(SBMLDocument, sbmlns)
        self.sbml_model = wrap_libsbml(self.sbml_document.createModel)

        self.per_second_id = 'per_second'
        self.per_second = wrap_libsbml(self.sbml_model.createUnitDefinition)
        wrap_libsbml(self.per_second.setIdAttribute, self.per_second_id)
        add_sbml_unit(self.per_second, UNIT_KIND_SECOND, exponent=-1)

    def test_add_sbml_unit(self):
        per_second = wrap_libsbml(self.sbml_model.createUnitDefinition)
        wrap_libsbml(per_second.setIdAttribute, 'per_second')
        self.assertTrue(wrap_libsbml(per_second.hasRequiredAttributes))
        exp = -1
        default_scale=0
        default_multiplier=1.0
        unit = add_sbml_unit(per_second, UNIT_KIND_SECOND, exponent=exp)
        self.assertEqual(wrap_libsbml(unit.getExponent, returns_int=True), exp)
        self.assertEqual(wrap_libsbml(unit.getKind), UNIT_KIND_SECOND)
        self.assertEqual(wrap_libsbml(unit.getScale), default_scale)
        self.assertEqual(wrap_libsbml(unit.getMultiplier), default_multiplier)

        strange_unit = wrap_libsbml(self.sbml_model.createUnitDefinition)
        wrap_libsbml(strange_unit.setIdAttribute, 'strange_unit')
        self.assertTrue(wrap_libsbml(strange_unit.hasRequiredAttributes))
        exp=-4; scale=3; mult=1.23
        unit = add_sbml_unit(strange_unit, UNIT_KIND_MOLE,
            exponent=exp, scale=scale, multiplier=mult)
        self.assertEqual(wrap_libsbml(unit.getExponent, returns_int=True), exp)
        self.assertEqual(wrap_libsbml(unit.getKind), UNIT_KIND_MOLE)
        self.assertEqual(wrap_libsbml(unit.getScale), scale)
        self.assertEqual(wrap_libsbml(unit.getMultiplier), mult)

        with self.assertRaises(LibSBMLError) as context:
            unit = add_sbml_unit(strange_unit, -1)
        self.assertIn("LibSBML returned error code", str(context.exception))

    def test_create_sbml_parameter(self):
        id='id1'; name='name1'; value=13; constant=False
        parameter = create_sbml_parameter(self.sbml_model, id, name=name, value=value, constant=constant)
        self.assertTrue(wrap_libsbml(parameter.hasRequiredAttributes))
        self.assertEqual(wrap_libsbml(parameter.getIdAttribute), id)
        self.assertEqual(wrap_libsbml(parameter.getName), name)
        self.assertTrue(wrap_libsbml(parameter.isSetValue))
        self.assertFalse(wrap_libsbml(parameter.isSetUnits))
        self.assertEqual(wrap_libsbml(parameter.getValue), value)
        self.assertEqual(wrap_libsbml(parameter.getConstant), constant)

        # test defaults
        id = 'id2'
        parameter = create_sbml_parameter(self.sbml_model, id)
        self.assertEqual(wrap_libsbml(parameter.getIdAttribute), id)
        self.assertEqual(wrap_libsbml(parameter.getName), '')
        self.assertFalse(wrap_libsbml(parameter.isSetValue))
        self.assertEqual(wrap_libsbml(parameter.getConstant), True)

        # test units
        id = 'id3'
        parameter = create_sbml_parameter(self.sbml_model, id, units=self.per_second_id)
        self.assertTrue(wrap_libsbml(parameter.hasRequiredAttributes))
        self.assertTrue(wrap_libsbml(parameter.isSetUnits))
        self.assertEqual(wrap_libsbml(parameter.getUnits), self.per_second_id)

        # test Parameter id collision
        with self.assertRaises(ValueError) as context:
            parameter = create_sbml_parameter(self.sbml_model, id)
        self.assertIn("is already in use as a Parameter id", str(context.exception))
