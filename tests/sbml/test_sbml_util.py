""" Tests of SBML utils

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import capturer
import unittest
import os
import six
import warnings

# "from libsbml import *" generates "NameError: Unknown C global variable" in pytest,
# presumably from the SWIG wrapper: http://web.mit.edu/svn/src/swig-1.3.25/Lib/python/pyinit.swg
from libsbml import (LIBSBML_OPERATION_SUCCESS, SBMLDocument, OperationReturnValue_toString,
                     UnitDefinition, SBMLNamespaces, UNIT_KIND_SECOND, UNIT_KIND_MOLE, UNIT_KIND_AMPERE,
                     UNIT_KIND_AVOGADRO)

from wc_lang.sbml.util import LibSbmlError, LibSbmlInterface


class TestSbml(unittest.TestCase):

    def setUp(self):
        # create an SBMLDocument that uses version 2 of the 'Flux Balance Constraints' extension
        self.document = LibSbmlInterface.create_doc(packages={'fbc': 2})

    def test_SBML_wrap_libsbml(self):

        id = u'test_id'
        self.assertEqual(
            LibSbmlInterface.call_libsbml(self.document.setIdAttribute, id), LIBSBML_OPERATION_SUCCESS)
        self.assertEqual(
            str(LibSbmlInterface.call_libsbml(self.document.getIdAttribute)), str(id))

        model = LibSbmlInterface.call_libsbml(self.document.createModel)
        self.assertEqual(
            LibSbmlInterface.call_libsbml(model.setTimeUnits, 'second'), LIBSBML_OPERATION_SUCCESS)

        self.assertEqual(
            LibSbmlInterface.call_libsbml(model.setTimeUnits, 'second',
                                          debug=True, returns_int=False), LIBSBML_OPERATION_SUCCESS)

        self.assertEqual(
            LibSbmlInterface.call_libsbml(self.document.getNumErrors, returns_int=True), 0)

        with self.assertRaises(LibSbmlError) as context:
            LibSbmlInterface.call_libsbml(self.document.getNumErrors, 'no arg')
        self.assertIn('Error', str(context.exception))
        self.assertIn("in libSBML method call", str(context.exception))

        with self.assertRaises(LibSbmlError) as context:
            LibSbmlInterface.call_libsbml(self.document.setIdAttribute, '..')
        self.assertIn('LibSBML returned error code', str(context.exception))
        self.assertIn("when executing", str(context.exception))

        with self.assertRaises(LibSbmlError) as context:
            LibSbmlInterface.call_libsbml(self.document.getAnnotation)
        self.assertIn('libSBML returned None when executing', str(context.exception))

    def test_returns_int(self):
        sbml_model = LibSbmlInterface.call_libsbml(self.document.createModel)
        avogadro_unit_def = LibSbmlInterface.call_libsbml(sbml_model.createUnitDefinition)
        avogadro_unit = LibSbmlInterface.add_unit(avogadro_unit_def, UNIT_KIND_AVOGADRO)

        # All 4 cases:
        # {returns_int == {T, F}} x {returns_int should be set in {T, F}}
        # 1. returns_int == True && returns_int should NOT be set
        # inconsequential because returns_int is not used
        self.assertEqual(
            LibSbmlInterface.call_libsbml(self.document.getModel, returns_int=True), sbml_model)

        # 2. returns_int == False && returns_int should NOT be set
        # correct usage
        # explicitly set returns_int = False
        self.assertEqual(
            LibSbmlInterface.call_libsbml(self.document.getModel, returns_int=False), sbml_model)
        # use default returns_int = False
        self.assertEqual(
            LibSbmlInterface.call_libsbml(self.document.getModel), sbml_model)

        # 3. returns_int == True && returns_int should be set
        # correct usage
        self.assertEqual(
            LibSbmlInterface.call_libsbml(avogadro_unit.getKind, returns_int=True, debug=True), UNIT_KIND_AVOGADRO)

        # 4. returns_int == False && returns_int should be set
        # 4a. libsbml call returns int == LIBSBML_OPERATION_SUCCESS
        # fails to raise warning because UNIT_KIND_AMPERE == LIBSBML_OPERATION_SUCCESS == 0
        # this situation cannot be avoided by wrap_libsbm, but causes no harm
        ampere_unit_def = LibSbmlInterface.call_libsbml(sbml_model.createUnitDefinition)
        ampere_unit = LibSbmlInterface.add_unit(ampere_unit_def, UNIT_KIND_AMPERE)
        self.assertEqual(
            LibSbmlInterface.call_libsbml(ampere_unit.getKind, returns_int=False), UNIT_KIND_AMPERE)

        # 4b. libsbml call returns int error code, NOT known by OperationReturnValue_toString()
        # raises warning because UNIT_KIND_AVOGADRO == 1 and LIBSBML_OPERATION_SUCCESS == 0
        with warnings.catch_warnings(record=True) as w:
            self.assertEqual(
                LibSbmlInterface.call_libsbml(avogadro_unit.getKind, returns_int=False), UNIT_KIND_AVOGADRO)
            self.assertEqual(len(w), 1)
            self.assertIn("Perhaps an integer value is being returned; ", str(w[-1].message))

        # 4c. libsbml call returns int error code, KNOWN by OperationReturnValue_toString()
        # raises warning in case the libsbml call returns an int value
        strange_unit = LibSbmlInterface.call_libsbml(sbml_model.createUnitDefinition)
        unit = LibSbmlInterface.call_libsbml(strange_unit.createUnit)
        with self.assertRaises(LibSbmlError) as context:
            unit_kind = -1
            LibSbmlInterface.call_libsbml(unit.setKind, unit_kind)
        self.assertIn("WARNING: if this libSBML call returns an int value, then this error may be incorrect",
                      str(context.exception))

    def test_init_sbml_model(self):
        sbml_model = LibSbmlInterface.init_model(None, self.document, packages={'fbc': 2})

        # check the SBML document
        self.assertEqual(LibSbmlInterface.call_libsbml(self.document.checkConsistency, returns_int=True), 0)
        self.assertTrue(LibSbmlInterface.is_doc_compatible(self.document))

        # check mmol_per_gDW_per_hr
        mmol_per_gDW_per_hr = LibSbmlInterface.call_libsbml(sbml_model.getUnitDefinition, 'mmol_per_gDW_per_hr')
        printed_mmol_per_gDW_per_hr = LibSbmlInterface.call_libsbml(UnitDefinition.printUnits, mmol_per_gDW_per_hr)
        compact_mmol_per_gDW_per_hr = LibSbmlInterface.call_libsbml(UnitDefinition.printUnits, mmol_per_gDW_per_hr, True)
        self.assertIn('(0.001 mole)^1', compact_mmol_per_gDW_per_hr)
        self.assertIn('(3600 second)^-1', compact_mmol_per_gDW_per_hr)
        self.assertIn('(1 gram)^-1', compact_mmol_per_gDW_per_hr)

    def test_SBML_fbc(self):

        # create an SBMLDocument that uses version 2 of the 'Flux Balance Constraints' extension
        document = LibSbmlInterface.create_doc(packages={'fbc': 2})

        id = 'x'
        self.assertEqual(
            LibSbmlInterface.call_libsbml(document.setIdAttribute, id), LIBSBML_OPERATION_SUCCESS)


class TestLibsbmlInterface(unittest.TestCase):

    def setUp(self):
        sbmlns = LibSbmlInterface.call_libsbml(SBMLNamespaces, 3, 1, "fbc", 1)
        self.sbml_document = LibSbmlInterface.call_libsbml(SBMLDocument, sbmlns)
        self.sbml_model = LibSbmlInterface.call_libsbml(self.sbml_document.createModel)

        self.per_second_id = 'per_second'
        self.per_second = LibSbmlInterface.call_libsbml(self.sbml_model.createUnitDefinition)
        LibSbmlInterface.call_libsbml(self.per_second.setIdAttribute, self.per_second_id)
        LibSbmlInterface.add_unit(self.per_second, UNIT_KIND_SECOND, exponent=-1)

    def test_add_sbml_unit(self):
        per_second = LibSbmlInterface.call_libsbml(self.sbml_model.createUnitDefinition)
        LibSbmlInterface.call_libsbml(per_second.setIdAttribute, 'per_second')
        self.assertTrue(LibSbmlInterface.call_libsbml(per_second.hasRequiredAttributes))
        exp = -1
        default_scale = 0
        default_multiplier = 1.0
        unit = LibSbmlInterface.add_unit(per_second, UNIT_KIND_SECOND, exponent=exp)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getExponent, returns_int=True), exp)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getKind, returns_int=True), UNIT_KIND_SECOND)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getScale, returns_int=True), default_scale)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getMultiplier), default_multiplier)

        strange_unit = LibSbmlInterface.call_libsbml(self.sbml_model.createUnitDefinition)
        LibSbmlInterface.call_libsbml(strange_unit.setIdAttribute, 'strange_unit')
        self.assertTrue(LibSbmlInterface.call_libsbml(strange_unit.hasRequiredAttributes))
        exp = -4
        scale = 3
        mult = 1.23
        unit = LibSbmlInterface.add_unit(strange_unit, UNIT_KIND_MOLE,
                                         exponent=exp, scale=scale, multiplier=mult)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getExponent, returns_int=True), exp)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getKind, returns_int=True), UNIT_KIND_MOLE)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getScale, returns_int=True), scale)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getMultiplier), mult)

        with self.assertRaises(LibSbmlError) as context:
            unit = LibSbmlInterface.add_unit(strange_unit, -1)
        self.assertIn("LibSBML returned error code", str(context.exception))

    def test_create_sbml_parameter(self):
        id = 'id1'
        name = 'name1'
        value = 13.0
        constant = False
        parameter = LibSbmlInterface.create_parameter(self.sbml_model, id, value, 'dimensionless_ud', name=name,
                                                      constant=constant)
        self.assertTrue(LibSbmlInterface.call_libsbml(parameter.hasRequiredAttributes))
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getIdAttribute), id)
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getName), name)
        self.assertTrue(LibSbmlInterface.call_libsbml(parameter.isSetValue))
        self.assertTrue(LibSbmlInterface.call_libsbml(parameter.isSetUnits))
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getValue), value)
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getConstant), constant)

        # test defaults
        id = 'id2'
        parameter = LibSbmlInterface.create_parameter(self.sbml_model, id, value, 'dimensionless_ud')
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getIdAttribute), id)
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getName), '')
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getConstant), True)

        # test units
        id = 'id3'
        parameter = LibSbmlInterface.create_parameter(self.sbml_model, id, value, units=self.per_second_id)
        self.assertTrue(LibSbmlInterface.call_libsbml(parameter.hasRequiredAttributes))
        self.assertTrue(LibSbmlInterface.call_libsbml(parameter.isSetUnits))
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getUnits), self.per_second_id)

        # test Parameter id collision
        with self.assertRaises(ValueError) as context:
            parameter = LibSbmlInterface.create_parameter(self.sbml_model, id, value, units=self.per_second_id)
        self.assertIn("is already in use as a Parameter id", str(context.exception))


class TestDebug(unittest.TestCase):

    def setUp(self):
        # create an SBMLDocument that uses version 2 of the 'Flux Balance Constraints' extension
        self.document = LibSbmlInterface.create_doc(packages={'fbc': 2})

    def test_SBML_wrap_libsbml_with_debug(self):
        with capturer.CaptureOutput() as capture_output:
            LibSbmlInterface.call_libsbml(self.document.setIdAttribute, 'id', debug=True)
            self.assertRegex(capture_output.get_text(), 'libSBML call:')
            self.assertRegex(capture_output.get_text(), 'libSBML returns:')

        sbml_model = LibSbmlInterface.call_libsbml(self.document.createModel)
        unit_def = LibSbmlInterface.call_libsbml(sbml_model.createUnitDefinition)
        unit = LibSbmlInterface.add_unit(unit_def, UNIT_KIND_AVOGADRO)
        with capturer.CaptureOutput() as capture_output:
            LibSbmlInterface.call_libsbml(unit.getKind, returns_int=False, debug=True)
            self.assertRegex(capture_output.get_text(), 'libSBML returns:')

        with capturer.CaptureOutput() as capture_output:
            LibSbmlInterface.call_libsbml(self.document.getIdAttribute, debug=True)
            self.assertRegex(capture_output.get_text(), 'libSBML returns:')
