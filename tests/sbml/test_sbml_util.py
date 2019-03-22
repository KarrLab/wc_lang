""" Tests of SBML utils

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2019-03-21
:Copyright: 2017-2019, Karr Lab
:License: MIT
"""

import capturer
import libsbml
import mock
import unittest
import warnings


# "from libsbml import *" generates "NameError: Unknown C global variable" in pytest,
# presumably from the SWIG wrapper: http://web.mit.edu/svn/src/swig-1.3.25/Lib/python/pyinit.swg
from libsbml import (LIBSBML_OPERATION_SUCCESS, SBMLDocument, OperationReturnValue_toString,
                     UnitDefinition, SBMLNamespaces, UNIT_KIND_SECOND, UNIT_KIND_MOLE, UNIT_KIND_AMPERE,
                     UNIT_KIND_AVOGADRO)

from wc_lang.core import Model, Species, ObservableExpression, WcLangWarning
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
        avogadro_unit = LibSbmlInterface.create_base_unit('avogadro_unit_def', avogadro_unit_def, 'avogadro')

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
        ampere_unit = LibSbmlInterface.create_base_unit('ampere_unit_def', ampere_unit_def, 'ampere')
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
        LibSbmlInterface.verify_doc_is_compatible(self.document)

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
        sbmlns = LibSbmlInterface.call_libsbml(SBMLNamespaces, 3, 2, "fbc", 2)
        self.sbml_document = LibSbmlInterface.call_libsbml(SBMLDocument, sbmlns)
        LibSbmlInterface.call_libsbml(self.sbml_document.setPackageRequired, 'fbc', False)
        self.sbml_model = LibSbmlInterface.call_libsbml(self.sbml_document.createModel)
        plugin = LibSbmlInterface.call_libsbml(self.sbml_model.getPlugin, 'fbc')
        LibSbmlInterface.call_libsbml(plugin.setStrict, True)

    def test_verify_doc_is_consistent(self):
        self.setUp()
        self.sbml_model.setTimeUnits('second')
        LibSbmlInterface.verify_doc_is_consistent(self.sbml_document)

        self.setUp()
        self.sbml_model.setTimeUnits('metre')
        with self.assertRaisesRegex(LibSbmlError, 'inconsistent'):
            LibSbmlInterface.verify_doc_is_consistent(self.sbml_document, strict_units=True)

        self.setUp()
        self.sbml_model.setTimeUnits('hour')
        with self.assertRaisesRegex(LibSbmlError, 'inconsistent'):
            LibSbmlInterface.verify_doc_is_consistent(self.sbml_document, strict_units=False)

    def test_verify_doc(self):
        LibSbmlInterface.verify_doc(self.sbml_document)

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
        unit = LibSbmlInterface.create_base_unit('per_second', per_second, 'second', exponent=exp)
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
        unit = LibSbmlInterface.create_base_unit('strange_unit', strange_unit, 'mole',
                                                 exponent=exp, scale=scale, multiplier=mult)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getExponent, returns_int=True), exp)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getKind, returns_int=True), UNIT_KIND_MOLE)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getScale, returns_int=True), scale)
        self.assertEqual(LibSbmlInterface.call_libsbml(unit.getMultiplier), mult)

        with self.assertRaisesRegex(AttributeError, 'no attribute'):
            LibSbmlInterface.create_base_unit('strange_unit', strange_unit, 'NOT_A_UNIT')

    def test_gen_unit_id(self):
        self.assertEqual(LibSbmlInterface.gen_unit_id(unit_registry.parse_units('meter / second')), 'unit_meter_per_second')
        self.assertEqual(LibSbmlInterface.gen_unit_id(unit_registry.parse_units('second')), 'second')

        with self.assertRaisesRegex(ValueError, 'Cannot generate SBML id'):
            LibSbmlInterface.gen_unit_id(None)

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
            unit_registry.parse_units('gDCW / s'),
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
            'unit_gDCW_per_second': units[6],
        }
        parsed_units = LibSbmlInterface.parse_units(self.sbml_model)
        self.assertEqual(parsed_units.keys(), exp_units.keys())
        for key in exp_units:
            self.assertTrue(are_units_equivalent(parsed_units[key], exp_units[key]),
                            '{} != {}'.format(str(parsed_units[key]), str(exp_units[key])))

    def test_parse_units_error(self):
        with self.assertRaisesRegex(LibSbmlError, 'units must be set'):
            LibSbmlInterface.parse_units(self.sbml_model)

    def test_create_parameter(self):
        self.per_second_id = 'unit_1_per_second'
        self.per_second = LibSbmlInterface.call_libsbml(self.sbml_model.createUnitDefinition)
        LibSbmlInterface.call_libsbml(self.per_second.setIdAttribute, self.per_second_id)
        LibSbmlInterface.create_base_unit('unit_1_per_second', self.per_second, 'second', exponent=-1)

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
        self.assertEqual(name, '')
        self.assertEqual(value, 1.5)
        self.assertEqual(units, unit_registry.parse_expression('s^-1'))

    def test_set_get_math(self):
        model = Model(version='1.2.3', wc_lang_version='4.5.6')

        spec_1_c = Species(id='spec_1[c]')
        spec_2_c = Species(id='spec_2[c]')
        model_objs = {Species: {spec_1_c.id: spec_1_c, spec_2_c.id: spec_2_c}}

        expression, error = ObservableExpression.deserialize('spec_1[c] + spec_2[c]', model_objs)
        assert error is None, str(error)

        sbml_doc = LibSbmlInterface.create_doc()
        sbml_model = LibSbmlInterface.init_model(model, sbml_doc)
        rule = LibSbmlInterface.call_libsbml(sbml_model.createAssignmentRule)
        LibSbmlInterface.set_math(rule.setMath, expression)

        self.assertEqual(libsbml.formulaToL3String(LibSbmlInterface.call_libsbml(rule.getMath)),
                         'Species__spec_1__RB__c__LB__ + Species__spec_2__RB__c__LB__')

        expression_2 = LibSbmlInterface.get_math(rule.getMath, ObservableExpression, model_objs)
        self.assertTrue(expression_2.is_equal(expression))
        self.assertEqual(expression_2._parsed_expression._obj_model_tokens,
                         expression._parsed_expression._obj_model_tokens)

    def test_export_import_annotations(self):
        model = Model()
        sbml_doc = LibSbmlInterface.create_doc()
        sbml_model = LibSbmlInterface.init_model(model, sbml_doc)
        LibSbmlInterface.set_annotations(model, {'version': 'version', 'wc_lang_version': 'wc_lang_version'}, sbml_model)

        model_2 = Model()
        LibSbmlInterface.get_annotations(model_2, {'version': 'version', 'wc_lang_version': 'wc_lang_version'}, sbml_model)
        self.assertEqual(model_2.version, model.version)
        self.assertEqual(model_2.wc_lang_version, model.wc_lang_version)

        model_3 = Model()
        LibSbmlInterface.get_annotations(model_3, {'version': [('version',)], 'wc_lang_version': [('wc_lang_version', )]}, sbml_model)
        self.assertEqual(model_3.version, model.version)
        self.assertEqual(model_3.wc_lang_version, model.wc_lang_version)

    def test_gen_nested_attr_paths(self):
        self.assertEqual(LibSbmlInterface.gen_nested_attr_paths(['version', 'wc_lang_version']),
                         {'version': [('version',)], 'wc_lang_version': [('wc_lang_version', )]})

        self.assertEqual(LibSbmlInterface.gen_nested_attr_paths(['model.version', 'model.wc_lang_version']),
                         {'model.version': [('model', ), ('version',)], 'model.wc_lang_version': [('model', ), ('wc_lang_version', )]})

    def test_str_to_xml_node(self):
        node = LibSbmlInterface.str_to_xml_node('')
        self.assertEqual(node.getNumChildren(), 0)

        node = LibSbmlInterface.str_to_xml_node('test')
        self.assertEqual(node.getNumChildren(), 1)
        self.assertEqual(node.getChild(0).getCharacters(), 'test')

        node = LibSbmlInterface.str_to_xml_node('&lt;test')
        self.assertEqual(node.getNumChildren(), 1)
        self.assertEqual(node.getChild(0).toXMLString(), '&lt;test')
        self.assertEqual(node.getChild(0).getCharacters(), '<test')

    def test_str_to_from_xml_node(self):
        for text in ['test', 'line1\nline2', 'line1<br/>line2', '<p>line1</p>\n<p>line2</p>']:
            node = LibSbmlInterface.str_to_xml_node(text)
            self.assertEqual(LibSbmlInterface.str_from_xml_node(node), text)

    def test_str_to_xml_node_error(self):
        with self.assertRaisesRegex(LibSbmlError, 'libSBML returned None'):
            with mock.patch('libsbml.XMLNode.convertStringToXMLNode', return_value=None):
                LibSbmlInterface.str_to_xml_node('test')

    def test_export_import_comments(self):
        for comments in ['', 'My comments', 'My\ncomments', 'My<br/>comments', '<p>My</p>\n<p>comments</p>']:
            sbml_doc = LibSbmlInterface.create_doc()
            sbml_model = LibSbmlInterface.create_model(sbml_doc)

            model = Model(comments=comments)
            LibSbmlInterface.set_commments(model, sbml_model)

            model_2 = Model()
            LibSbmlInterface.get_commments(model_2, sbml_model)

            self.assertEqual(model_2.comments, comments)


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
        unit = LibSbmlInterface.create_base_unit('avogadro_unit_def', unit_def, 'avogadro')
        with self.assertWarnsRegex(WcLangWarning, 'unknown error code'):
            with capturer.CaptureOutput() as capture_output:
                LibSbmlInterface.call_libsbml(unit.getKind, returns_int=False, debug=True)
                self.assertRegex(capture_output.get_text(), 'libSBML returns:')

        with capturer.CaptureOutput() as capture_output:
            LibSbmlInterface.call_libsbml(self.document.getIdAttribute, debug=True)
            self.assertRegex(capture_output.get_text(), 'libSBML returns:')
