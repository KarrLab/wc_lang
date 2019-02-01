""" Tests of SBML utils

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import capturer
import mock
import os
import six
import unittest
import warnings


# "from libsbml import *" generates "NameError: Unknown C global variable" in pytest,
# presumably from the SWIG wrapper: http://web.mit.edu/svn/src/swig-1.3.25/Lib/python/pyinit.swg
from libsbml import (LIBSBML_OPERATION_SUCCESS, SBMLDocument, OperationReturnValue_toString,
                     UnitDefinition, SBMLNamespaces, UNIT_KIND_SECOND, UNIT_KIND_MOLE, UNIT_KIND_AMPERE,
                     UNIT_KIND_AVOGADRO)

from wc_lang.core import Model
from wc_lang.sbml.util import LibSbmlError, LibSbmlInterface
from wc_utils.util.units import unit_registry, are_units_equivalent


class TestSbml(unittest.TestCase):

    def setUp(self):
        # create an SBMLDocument that uses version 2 of the 'Flux Balance Constraints' extension
        self.document = LibSbmlInterface.create_doc(packages={'fbc': 2})

    def test_LibSbmlError(self):
        error = LibSbmlError('test')
        self.assertEqual(str(error), 'test')

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
        avogadro_unit = LibSbmlInterface.create_base_unit(avogadro_unit_def, 'avogadro')

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
        ampere_unit = LibSbmlInterface.create_base_unit(ampere_unit_def, 'ampere')
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

    def test_init_model(self):
        model = Model()
        model.parameters.create(units=unit_registry.parse_units('s^-1'))
        model.parameters.create(units=unit_registry.parse_units('g / l'))
        model.parameters.create(units=unit_registry.parse_units('m'))
        sbml_model = LibSbmlInterface.init_model(model, self.document, packages={'fbc': 2})

        # check the SBML document
        return_val = LibSbmlInterface.call_libsbml(self.document.checkConsistency, returns_int=True)
        self.assertEqual(return_val, 0)
        self.assertTrue(LibSbmlInterface.is_doc_compatible(self.document))

        # check seconds unit
        unit_def = LibSbmlInterface.call_libsbml(sbml_model.getUnitDefinition, 'unit_1_per_second')
        units = LibSbmlInterface.call_libsbml(UnitDefinition.printUnits, unit_def, True)
        self.assertEqual(units, '(1 second)^-1')

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

    def test_gen_base_unit_id(self):
        self.assertEqual(LibSbmlInterface.gen_base_unit_id(unit_registry.parse_units('mole')), 'mole')
        self.assertEqual(LibSbmlInterface.gen_base_unit_id(unit_registry.parse_units('meter')), 'metre')
        self.assertEqual(LibSbmlInterface.gen_base_unit_id(unit_registry.parse_units('liter')), 'litre')

    def normalize_unit_kind(self):
        self.assertEqual(LibSbmlInterface.normalize_unit_kind(unit_registry.parse_units('meter')), 'metre')
        self.assertEqual(LibSbmlInterface.normalize_unit_kind(unit_registry.parse_units('liter')), 'litre')
        self.assertEqual(LibSbmlInterface.normalize_unit_kind(unit_registry.parse_units('mole')), 'mole')
        self.assertEqual(LibSbmlInterface.normalize_unit_kind(unit_registry.parse_units('item')), 'item')

    def test_create_unit(self):
        per_second = LibSbmlInterface.call_libsbml(self.sbml_model.createUnitDefinition)
        LibSbmlInterface.call_libsbml(per_second.setIdAttribute, 'per_second')
        self.assertTrue(LibSbmlInterface.call_libsbml(per_second.hasRequiredAttributes))
        exp = -1
        default_scale = 0
        default_multiplier = 1.0
        unit = LibSbmlInterface.create_base_unit(per_second, 'second', exponent=exp)
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
        unit = LibSbmlInterface.create_base_unit(strange_unit, 'mole',
                                                 exponent=exp, scale=scale, multiplier=mult)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getExponent, returns_int=True), exp)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getKind, returns_int=True), UNIT_KIND_MOLE)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getScale, returns_int=True), scale)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getMultiplier), mult)

        with self.assertRaisesRegex(AttributeError, 'no attribute'):
            LibSbmlInterface.create_base_unit(strange_unit, 'NOT_A_UNIT')

    def test_parse_units(self):
        self.sbml_model.setTimeUnits('second')
        self.sbml_model.setSubstanceUnits('item')
        self.sbml_model.setExtentUnits('item')
        self.sbml_model.setVolumeUnits('litre')
        units = [
            unit_registry.parse_units('g'),
            unit_registry.parse_units('l'),
            unit_registry.parse_units('g/l'),
            unit_registry.parse_units('g/s'),
            unit_registry.parse_units('g/s**2'),
            unit_registry.parse_units('g/(l * ms)^2'),
        ]
        for unit in units:
            LibSbmlInterface.create_unit(unit, self.sbml_model)

        exp_units = {
            'second': unit_registry.parse_units('second'),
            'item': unit_registry.parse_units('molecule'),
            'litre': unit_registry.parse_units('liter'),
            'unit_gram_per_liter': units[2],
            'unit_gram_per_second': units[3],
            'unit_gram_per_second_pow_2': units[4],
            'unit_gram_per_liter_pow_2_per_millisecond_pow_2': units[5],
        }
        parsed_units = LibSbmlInterface.parse_units(self.sbml_model)
        self.assertEqual(parsed_units.keys(), exp_units.keys())
        for key in exp_units:
            self.assertTrue(are_units_equivalent(parsed_units[key], exp_units[key]))

    def test_parse_units_error(self):
        with self.assertRaisesRegex(LibSbmlError, 'units must be set'):
            LibSbmlInterface.parse_units(self.sbml_model)

    def test_create_parameter(self):
        self.per_second_id = 'unit_1_per_second'
        self.per_second = LibSbmlInterface.call_libsbml(self.sbml_model.createUnitDefinition)
        LibSbmlInterface.call_libsbml(self.per_second.setIdAttribute, self.per_second_id)
        LibSbmlInterface.create_base_unit(self.per_second, 'second', exponent=-1)

        id = 'id1'
        name = 'name1'
        value = 13.0
        constant = False
        parameter = LibSbmlInterface.create_parameter(self.sbml_model, id, value,
                                                      unit_registry.parse_units('dimensionless'), name=name,
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
        parameter = LibSbmlInterface.create_parameter(self.sbml_model, id, value,
                                                      unit_registry.parse_units('dimensionless'))
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getIdAttribute), id)
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getName), '')
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getConstant), True)

        # test units
        id = 'id3'
        parameter = LibSbmlInterface.create_parameter(self.sbml_model, id, value,
                                                      unit_registry.parse_units('s^-1'))
        self.assertTrue(LibSbmlInterface.call_libsbml(parameter.hasRequiredAttributes))
        self.assertTrue(LibSbmlInterface.call_libsbml(parameter.isSetUnits))
        self.assertEqual(LibSbmlInterface.call_libsbml(parameter.getUnits), self.per_second_id)

    def test_parse_parameter(self):
        units = unit_registry.parse_units('s^-1')
        LibSbmlInterface.create_unit(units, self.sbml_model)
        parameter = LibSbmlInterface.create_parameter(self.sbml_model, 'parameter_1', 1.5, units)

        id, name, value, units = LibSbmlInterface.parse_parameter(parameter)
        self.assertEqual(id, 'parameter_1')
        self.assertEqual(name, None)
        self.assertEqual(value, 1.5)
        self.assertEqual(units, unit_registry.parse_expression('s^-1'))

    def test_parse_parameter_error(self):
        units = unit_registry.parse_units('s^-1')
        LibSbmlInterface.create_unit(units, self.sbml_model)

        parameter = LibSbmlInterface.create_parameter(self.sbml_model, 'parameter_1', 1.5, units,
                                                      constant=False)
        with self.assertRaisesRegex(LibSbmlError, 'must be constant'):
            LibSbmlInterface.parse_parameter(parameter)

        units = unit_registry.parse_units('g^-1')
        parameter = LibSbmlInterface.create_parameter(self.sbml_model, 'parameter_1', 1.5, units,
                                                      constant=False)
        with self.assertRaisesRegex(LibSbmlError, 'must be defined'):
            LibSbmlInterface.parse_parameter(parameter)

    def test_str_to_xml_node(self):
        node = LibSbmlInterface.str_to_xml_node('test')
        self.assertEqual(node.getNumChildren(), 1)
        self.assertEqual(node.getChild(0).getCharacters(), 'test')

        node = LibSbmlInterface.str_to_xml_node('<test')
        self.assertEqual(node.getNumChildren(), 1)
        self.assertEqual(node.getChild(0).getCharacters(), '<test')

    def test_str_to_xml_node_error(self):
        with self.assertRaisesRegex(LibSbmlError, 'libSBML returned None'):
            with mock.patch('libsbml.XMLNode.convertStringToXMLNode', return_value=None):
                LibSbmlInterface.str_to_xml_node('test')


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
        unit = LibSbmlInterface.create_base_unit(unit_def, 'avogadro')
        with capturer.CaptureOutput() as capture_output:
            LibSbmlInterface.call_libsbml(unit.getKind, returns_int=False, debug=True)
            self.assertRegex(capture_output.get_text(), 'libSBML returns:')

        with capturer.CaptureOutput() as capture_output:
            LibSbmlInterface.call_libsbml(self.document.getIdAttribute, debug=True)
            self.assertRegex(capture_output.get_text(), 'libSBML returns:')
