""" Utilities for writing/reading a `wc_lang` model to/from SBML

* Higher level functions for creating, getting, and setting SBML objects
* Utilities for wrapping libSBML calls

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2019-03-21
:Copyright: 2017-2019, Karr Lab
:License: MIT
"""

from wc_utils.util.units import unit_registry, are_units_equivalent
import libsbml
import math
import obj_model
import obj_model.expression
import obj_model.units
import re
import scipy
import types
import warnings
import wc_lang.core


class SbmlModelMixin(object):
    """ Methods for exporting/import models to/from SBML """

    def gen_sbml_id(self):
        """ Generate SBML id from id

        Returns:
            :obj:`str`: SBML id
        """
        return self.__class__.__name__ + '__' \
            + self.id \
            .replace('[', '__RB__') \
            .replace(']', '__LB__') \
            .replace('-', '__DS__')

    @classmethod
    def parse_sbml_id(cls, sbml_id):
        """ Parse id from SBML id

        Args:
            sbml_id (:obj:`str`): SBML id

        Returns:
            :obj:`str`: id
        """
        return sbml_id.partition(cls.__name__ + '__')[2] \
            .replace('__RB__', '[') \
            .replace('__LB__', ']') \
            .replace('__DS__', '-')

    def export_to_sbml(self, sbml_model):
        """ Add object to SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.SBase`: SBML object
        """
        pass  # pragma: no cover

    def export_relations_to_sbml(self, sbml_model, sbml):
        """ Add relationships to/from object to SBML object.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml (:obj:`libsbml.SBase`): SBML object
        """
        pass  # pragma: no cover

    def import_from_sbml(self, sbml):
        """ Load object from SBML object

        Args:
            sbml (:obj:`libsbml.SBase`): SBML object
        """
        pass  # pragma: no cover

    def import_relations_from_sbml(self, sbml, objs):
        """ Load relations to/from object from SBML object

        Args:
            sbml (:obj:`libsbml.SBase`): SBML object
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        pass  # pragma: no cover


class SbmlAssignmentRuleMixin(SbmlModelMixin):
    """ Methods for exporting/import expression models to/from SBML """

    def export_to_sbml(self, sbml_model):
        """ Add expression to a SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.AssignmentRule`: SBML assignment rule
        """
        sbml = LibSbmlInterface.call_libsbml(sbml_model.createAssignmentRule)

        LibSbmlInterface.call_libsbml(sbml.setIdAttribute, self.gen_sbml_id())
        LibSbmlInterface.call_libsbml(sbml.setName, self.name)

        param_id = '__param__' + self.gen_sbml_id()
        LibSbmlInterface.create_parameter(sbml_model, param_id, None, self.units, constant=False)
        LibSbmlInterface.call_libsbml(sbml.setVariable, param_id)

        LibSbmlInterface.set_commments(self, sbml)

        return sbml

    def export_relations_to_sbml(self, sbml_model, sbml):
        """ Add relationships to/from object to SBML object.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml (:obj:`libsbml.SBase`): SBML object
        """
        LibSbmlInterface.set_math(sbml.setMath, self.expression)
        LibSbmlInterface.set_annotations(self, LibSbmlInterface.gen_nested_attr_paths(['identifiers']), sbml)

    def import_from_sbml(self, sbml):
        """ Load expression from SBML assignment rule

        Args:
            sbml (:obj:`libsbml.AssignmentRule`): SBML assignment rule
        """
        self.id = self.parse_sbml_id(LibSbmlInterface.call_libsbml(sbml.getIdAttribute))
        self.name = LibSbmlInterface.call_libsbml(sbml.getName)

        param_id = LibSbmlInterface.call_libsbml(sbml.getVariable)
        sbml_model = LibSbmlInterface.call_libsbml(sbml.getModel)
        _, _, _, self.units = LibSbmlInterface.parse_parameter(
            LibSbmlInterface.call_libsbml(sbml_model.getParameter, param_id))

        LibSbmlInterface.get_commments(self, sbml)

    def import_relations_from_sbml(self, sbml, objs):
        """ Load relationships to/from expression from SBML assignment rule

        Args:
            sbml (:obj:`libsbml.AssignmentRule`): SBML assignment rule
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        self.expression = LibSbmlInterface.get_math(sbml.getMath, self.Meta.expression_term_model, objs)
        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(['identifiers']), sbml, objs)


class LibSbmlError(Exception):
    ''' Exception raised when libSBML returns an error '''


class LibSbmlInterface(object):
    ''' Methods for compactly using libSBML to create SBML objects.

    The libSBML method calls provide narrow interfaces, typically exchanging one
    value per call, which creates verbose code. The methods below aggregate multiple
    libSBML method calls to enable more compact usage.
    '''

    XML_NAMESPACE = 'https://www.wholecell.org/ns/wc_lang'

    @classmethod
    def create_doc(cls, level=3, version=2, packages=None):
        """ Create an SBMLDocument that, optionally, uses package(s).

        Args:
            level (:obj:`int`, optional): SBML level number
            version (:obj:`int`, optional): SBML version number
            packages (:obj:`dict` that maps :obj:`str` to :obj:`int`, optional): dictionary of required packages
                that maps package identifiers to package numbers

        Returns:
            :obj:`libsbml.SBMLDocument`: SBML document
        """
        packages = packages or {}

        # create name spaces for SBML document
        namespaces = [level, version]
        for package_id, package_version in packages.items():
            namespaces.append(package_id)
            namespaces.append(package_version)
        sbml_ns = cls.call_libsbml(libsbml.SBMLNamespaces, *namespaces)
        cls.call_libsbml(sbml_ns.addNamespace, cls.XML_NAMESPACE, 'wcLang')

        # create SBML document
        sbml_doc = cls.call_libsbml(libsbml.SBMLDocument, sbml_ns)

        # set package requirements
        for package in packages:
            cls.call_libsbml(sbml_doc.setPackageRequired, package, False)

        # return SBML document
        return sbml_doc

    @classmethod
    def verify_doc(cls, sbml_doc, level=3, version=2, strict_units=True):
        """ Verify that an SBML document is valid, and raise an exception if the document is not valid

        * SBML-compatible
        * Valid SBML
        * Consistent

        Args:
            sbml_doc (:obj:`libsbml.SBMLDocument`): SBML document
            level (:obj:`int`, optional): SBML level number
            version (:obj:`int`, optional): SBML version number
            strict_units (:obj:`bool`, optional): if :obj:`True`, strictly verify that the units are
                consistent
        """
        cls.verify_doc_is_compatible(sbml_doc, level=level, version=version)
        cls.verify_doc_is_valid_sbml(sbml_doc)
        cls.verify_doc_is_consistent(sbml_doc, strict_units=strict_units)

    @classmethod
    def verify_doc_is_compatible(cls, sbml_doc, level=3, version=2):
        """ Check the compatibility of an SBML document with a specific level and version, 
        and raise an exception if the document is not verify_doc_is_compatible

        Args:
            sbml_doc (:obj:`libsbml.SBMLDocument`): SBML document
            level (:obj:`int`, optional): SBML level number
            version (:obj:`int`, optional): SBML version number
        """
        # SBML compatibility method for the version being used
        method = getattr(sbml_doc, 'checkL{}v{}Compatibility'.format(level, version))
        cls.call_libsbml(method, returns_int=True)
        cls.raise_if_error(sbml_doc, 'Document is incompatible')

    @classmethod
    def verify_doc_is_valid_sbml(cls, sbml_doc):
        """ Check that an SBML document is valid, 
        and raise an exception if the document is not valid

        Args:
            sbml_doc (:obj:`libsbml.SBMLDocument`): SBML document
        """
        cls.call_libsbml(sbml_doc.validateSBML, returns_int=True)
        cls.raise_if_error(sbml_doc, 'Document is invalid SBML')

    @classmethod
    def verify_doc_is_consistent(cls, sbml_doc, strict_units=True):
        """ Check that an SBML document is consistent, 
        and raise an exception if the document is not consistent

        Args:
            sbml_doc (:obj:`libsbml.SBMLDocument`): SBML document
            strict_units (:obj:`bool`, optional): if :obj:`True`, strictly verify that the units are
                consistent
        """
        cls.call_libsbml(sbml_doc.checkInternalConsistency, returns_int=True)
        cls.raise_if_error(sbml_doc, 'Document is internally inconsistent')

        if strict_units:
            method = sbml_doc.checkConsistencyWithStrictUnits
        else:
            method = sbml_doc.checkConsistency
        cls.call_libsbml(method, returns_int=True)
        cls.raise_if_error(sbml_doc, 'Document is inconsistent')

    @classmethod
    def create_model(cls, sbml_doc):
        """ Create a SBML model

        Args:
            sbml_doc (:obj:`libsbml.SBMLDocument`): SBML document

        Returns:
            :obj:`libsbml.Model`: SBML model
        """
        return cls.call_libsbml(sbml_doc.createModel)

    @classmethod
    def init_model(cls, model, sbml_doc, packages=None):
        """ Create and initialize an SMBL model.

        Args:
            model (:obj:`wc_lang.core.Model`): model
            sbml_doc (:obj:`libsbml.SBMLDocument`): a `libsbml` SBMLDocument
            packages (:obj:`dict` that maps :obj:`str` to :obj:`int`, optional): dictionary of required packages
                that maps package identifiers to package numbers

        Returns:
            :obj:`libsbml.Model`: the SBML model
        """

        # create model
        sbml_model = cls.create_model(sbml_doc)

        # enable plugins for packages
        packages = packages or {}
        for package_id in packages.keys():
            plugin = cls.call_libsbml(sbml_model.getPlugin, package_id)
            cls.call_libsbml(plugin.setStrict, True)

        # Set units
        cls.create_units(model, sbml_model)

        # return SBML model
        return sbml_model

    @classmethod
    def create_units(cls, model, sbml_model):
        """ Set time, extent, and substance units of SBML model

        Args:
            model (:obj:`wc_lang.core.Model`): model
            sbml_model (:obj:`libsbml.Model`): SBML model that encodes the model

        Returns:
            :obj:`dict`: dictionary that maps units to ids of SBML unit definitions
        """
        # Define units
        units = obj_model.units.get_obj_units(model)

        units_to_sbml = {}
        for unit in units:
            sbml_unit = cls.create_unit(unit, sbml_model)
            if sbml_unit:
                units_to_sbml[unit] = cls.call_libsbml(sbml_unit.getIdAttribute)

        # set top-level units (time, substance, extent, volume)
        cls.set_unit(sbml_model.setTimeUnits, model.time_units)

        assert len(wc_lang.core.Species.units.choices) == 1
        units = wc_lang.core.Species.units.choices[0]
        cls.set_unit(sbml_model.setExtentUnits, units)
        cls.set_unit(sbml_model.setSubstanceUnits, units)

        assert len(wc_lang.core.InitVolume.units.choices) == 1
        units = wc_lang.core.InitVolume.units.choices[0]
        cls.set_unit(sbml_model.setVolumeUnits, units)

        # return dictionary to SBML units
        return units_to_sbml

    @classmethod
    def create_unit(cls, unit, sbml_model):
        """ Add a unit definition to a SBML model

        Args:
            unit (:obj:`unit_registry.Unit`): unit
            sbml_model (:obj:`libsbml.Model`): SBML model that encodes the model

        Returns:
            :obj:`libsbml.UnitDefinition`: unit definition
        """
        id = cls.gen_unit_id(unit)
        if not id.startswith('unit_'):
            return None

        unit_def = cls.call_libsbml(sbml_model.createUnitDefinition)
        cls.call_libsbml(unit_def.setIdAttribute, id)

        unit_registry = unit._REGISTRY

        original_unit = unit

        quant = unit_registry.parse_expression(str(unit))

        mag, bases = quant.to_reduced_units().to_tuple()
        mole_bases = []
        has_molecule = False
        for base, exp in bases:
            if base == 'molecule':
                has_molecule = True
                base = 'mole'
            mole_bases.append((base, exp))
        mole_bases = tuple(mole_bases)
        if has_molecule:
            quant = unit_registry.Quantity.from_tuple((mag, mole_bases))

        magnitude, root_units = quant.to_root_units().to_tuple()
        scale = int(math.floor(math.log10(magnitude)))
        multiplier = magnitude / pow(10, scale)

        if not root_units:
            root_units = [('dimensionless', 1.), ]

        if are_units_equivalent(original_unit, unit_registry.parse_units('molecule / mole')):
            scale = 0
            multiplier = 1.
            root_units = [('dimensionless', 1.), ]

        elif are_units_equivalent(original_unit, unit_registry.parse_units('M')):
            scale = 0
            multiplier = 1.
            root_units = [
                ('mole', 1.),
                ('liter', -1.),
            ]

        for i_root_unit, (kind, exponent) in enumerate(root_units):
            if i_root_unit == 0:
                unit_scale = scale
                unit_multiplier = multiplier
            else:
                unit_scale = 0
                unit_multiplier = 1.

            cls.create_base_unit(id, unit_def, kind, exponent=exponent, scale=unit_scale, multiplier=unit_multiplier)

        return unit_def

    @classmethod
    def gen_unit_id(cls, unit):
        """ Generate an SBML unit id for a unit

        Args:
            unit (:obj:`unit_registry.Unit`): unit

        Returns:
            :obj:`str`: SBML id for unit

        Raises:
            :obj:`ValueError`: if an SBML id cannot be generated for the unit
        """
        if not isinstance(unit, unit_registry.Unit):
            raise ValueError('Cannot generate SBML id for `None` unit')

        id = 'unit_' + str(unit) \
            .replace(' / ', '_per_') \
            .replace(' * ', '_times_') \
            .replace(' ** ', '_pow_')

        if hasattr(libsbml, 'UNIT_KIND_' + cls.normalize_unit_kind(id[5:], to_sbml_base_units=False).upper()):
            id = cls.normalize_unit_kind(id[5:], to_sbml_base_units=False)

        return id

    @classmethod
    def create_base_unit(cls, unit_def_id, unit_def, kind, exponent=1, scale=0, multiplier=1.0):
        """ Add a unit to a SBML unit definition

        Each SBML unit has four attributes:

        * `kind`
        * `exponent`
        * `scale`
        * `multiplier`

        Args:
            unit_def_id (:obj:`str`): id of SBML unit definition
            unit_def (:obj:`libsbml.UnitDefinition`): SBML unit definition
            kind (:obj:`str`): unit kind
            exponent (:obj:`int`, optional): exponent of the unit
            scale (:obj:`int`, optional): scale of the unit
            multiplier (:obj:`float`, optional): multiplier of the unit

        Returns:
            :obj:`libsbml.Unit`: SBML unit
        """
        id = unit_def_id + '_' + kind
        kind = cls.normalize_unit_kind(kind)
        kind_val = getattr(libsbml, 'UNIT_KIND_' + kind.upper())

        unit = cls.call_libsbml(unit_def.createUnit)
        cls.call_libsbml(unit.setIdAttribute, id)
        cls.call_libsbml(unit.setKind, kind_val)
        cls.call_libsbml(unit.setExponent, exponent)
        cls.call_libsbml(unit.setScale, scale)
        cls.call_libsbml(unit.setMultiplier, multiplier)
        return unit

    @classmethod
    def normalize_unit_kind(cls, kind, to_sbml_base_units=True):
        """ Normalize unit kind for SBML

        Args:
            kind (:obj:`str`): unit kind
            to_sbml_base_units (:obj:`bool`, optional): if :obj:`True`, map to fundamental SBML units

        Returns:
            :obj:`str`: normalized unit kind
        """
        if kind == 'liter':
            kind = 'litre'
        elif kind == 'meter':
            kind = 'metre'
        elif kind == 'molecule':
            kind = 'mole'

        if to_sbml_base_units:
            if kind == 'gDCW':
                kind = 'gram'

        return kind

    @classmethod
    def parse_units(cls, sbml_model):
        """ Parse SBML units to Python

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`dict`: dictionary that maps ids of unit definitions to instance of `unit_registry.Unit`

        Raises:
            :obj:`LibSbmlError`: if units are not set or are inconsistent
        """
        units = {}

        for type in ['Time', 'Substance', 'Volume']:
            if not getattr(sbml_model, 'isSet' + type + 'Units')():
                raise LibSbmlError('{} units must be set'.format(type))

            sbml_unit_id = getattr(sbml_model, 'get' + type + 'Units')()
            if sbml_unit_id == 'mole':
                units[sbml_unit_id] = unit_registry.parse_expression('molecule')
            else:
                units[sbml_unit_id] = unit_registry.parse_expression(sbml_unit_id)

        assert cls.call_libsbml(sbml_model.getExtentUnits) == cls.call_libsbml(sbml_model.getSubstanceUnits), \
            LibSbmlError('Substance and extent units must be the same')
        assert not cls.call_libsbml(sbml_model.isSetLengthUnits), \
            LibSbmlError('Length units must be unset')
        assert not cls.call_libsbml(sbml_model.isSetAreaUnits), \
            LibSbmlError('Area units must be unset')
        for i_unit_def in range(cls.call_libsbml(sbml_model.getNumUnitDefinitions, returns_int=True)):
            sbml_unit_def = cls.call_libsbml(sbml_model.getUnitDefinition, i_unit_def)
            sbml_unit_def_id = cls.call_libsbml(sbml_unit_def.getIdAttribute)
            units[sbml_unit_def_id] = cls.parse_unit_id(sbml_unit_def_id)
        return units

    @classmethod
    def parse_unit_id(cls, sbml_unit_def_id):
        """ Parse a unit from an SBML unit definition.

        Args:
            sbml_unit_def_id (:obj:`str`): id of SBML unit definition

        Returns:
            :obj:`unit_registry.Unit`: units
        """
        if sbml_unit_def_id.startswith('unit_'):
            return unit_registry.parse_units(sbml_unit_def_id.partition('unit_')[2]
                                             .replace('_per_', ' / ')
                                             .replace('_times_', ' * ')
                                             .replace('_pow_', ' ** '))
        else:
            if sbml_unit_def_id == 'mole':
                return unit_registry.parse_units('molecule')
            return unit_registry.parse_units(sbml_unit_def_id)

    @classmethod
    def set_unit(cls, sbml_set_unit_func, unit):
        """ Set the units of an SBML object

        Args:
           sbml_set_unit_func (:obj:`callable`): function to set the units of an SBML object
           unit (:obj:`unit_registry.Unit`): unit
        """
        unit_id = cls.gen_unit_id(unit)
        cls.call_libsbml(sbml_set_unit_func, unit_id)

    @classmethod
    def get_unit(cls, sbml_get_unit_func):
        """ Get units from SBML unit definition id

        Args:
            sbml_get_unit_func (:obj:`callable`): function to get the units of an SBML object

        Returns:
            :obj:`unit_registry.Unit`: units
        """
        return cls.parse_unit_id(cls.call_libsbml(sbml_get_unit_func))

    @classmethod
    def create_parameter(cls, sbml_model, id, value, units, name=None, constant=True):
        """ Add a parameter to an SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            id (:obj:`str`): id
            value (:obj:`obj`): value
            units (:obj:`unit_registry.Unit`): units
            name (:obj:`str`, optional): name
            constant (:obj:`str`, optional): whether the parameter is a constant

        Returns:
            :obj:`libsbml.Parameter`: SBML parameter
        """
        sbml_parameter = cls.call_libsbml(sbml_model.createParameter)
        cls.call_libsbml(sbml_parameter.setIdAttribute, id)
        if name is not None:
            cls.call_libsbml(sbml_parameter.setName, name)
        if value is not None:
            cls.call_libsbml(sbml_parameter.setValue, value)
        cls.set_unit(sbml_parameter.setUnits, units)
        cls.call_libsbml(sbml_parameter.setConstant, constant)
        return sbml_parameter

    @classmethod
    def parse_parameter(cls, sbml_parameter):
        """ Parse a parameter from an SBML parameter.

        Args:
            sbml_parameter (:obj:`libsbml.Parameter`): SBML parameter

        Returns:
            :obj:`str`: id
            :obj:`str`: name
            :obj:`float`: value
            :obj:`unit_registry.Unit`: units
        """
        id = cls.call_libsbml(sbml_parameter.getIdAttribute)
        name = cls.call_libsbml(sbml_parameter.getName)
        value = cls.call_libsbml(sbml_parameter.getValue)
        units = cls.get_unit(sbml_parameter.getUnits)
        return (id, name, value, units)

    @classmethod
    def set_math(cls, set_math_func, expression, units_transform=None):
        """ Set the math of an SBML object

        Args:
            set_math_func (:obj:`callable`): function to set the math of an SBML object
            expression (:obj:`obj_model.expression.Expression`): expression
        """
        str_formula = expression._parsed_expression.get_str(cls._obj_model_token_to_str, with_units=True, number_units=' dimensionless')
        if units_transform:
            str_formula = units_transform.format(str_formula)
        sbml_formula = cls.call_libsbml(libsbml.parseL3Formula, str_formula)
        cls.call_libsbml(set_math_func, sbml_formula)

    @staticmethod
    def _obj_model_token_to_str(token):
        """ Get a string representation of a token that represents an instance of :obj:`Model`.

        Args:
            token (:obj:`obj_model.expression.ObjModelToken`): token that represents an instance of :obj:`Model`

        Returns:
            :obj:`str`: string representation of a token that represents an instance of :obj:`Model`.
        """
        if isinstance(token.model, wc_lang.core.Compartment):
            return '__mass__' + token.model.gen_sbml_id()
        elif isinstance(token.model, (wc_lang.core.Observable, wc_lang.core.Function)):
            return '__param__' + token.model.gen_sbml_id()
        else:
            return token.model.gen_sbml_id()

    @classmethod
    def get_math(cls, get_math_func, Expression, model_objs, units_transform=None):
        """ Get the math of an SBML object

        Args:
            get_math_func (:obj:`callable`): function to get the math of an SBML object
            Expression (:obj:`type`): type of expression
            model_objs (:obj:`dict`, optional): dictionary that maps classes of model objects to dictonaries
                that map ids of model objects to model objects

        Returns:
            :obj:`obj_model.expression.Expression`: expression
        """
        sbml_formula = cls.call_libsbml(get_math_func)
        str_formula = cls.call_libsbml(libsbml.formulaToL3String, sbml_formula)
        if units_transform:
            str_formula = units_transform(str_formula)
        str_formula = str_formula \
            .replace('__mass__Compartment__', 'Compartment__') \
            .replace('__param__Observable__', 'Observable__') \
            .replace('__param__Function__', 'Function__') \
            .replace(' dimensionless', '')
        str_formula = str_formula \
            .replace('__RB__', '[') \
            .replace('__LB__', ']') \
            .replace('__DS__', '-')
        str_formula = re.sub(r'[A-Za-z0-9]+__', '', str_formula)
        expression, error = Expression.deserialize(str_formula, model_objs)
        assert error is None, str(error)
        return expression

    @classmethod
    def set_annotations(cls, model_obj, nested_attr_paths, sbml_obj):
        """ Export annotations from a model object to an SBML object

        Args:
            model_obj (:obj:`obj_model.Model`): model object
            nested_attr_paths (:obj:`dict`): dictionary that maps names of attributes to paths to the attributes
            sbml_obj (:obj:`libsbml.SBase`): SBML object
        """
        cls.call_libsbml(sbml_obj.setAnnotation,
                         '<annotation><wcLang:annotation>'
                         + cls.gen_annotations(model_obj, nested_attr_paths, sbml_obj)
                         + '</wcLang:annotation></annotation>')

    @classmethod
    def gen_annotations(cls, model_obj, nested_attr_paths, sbml_obj):
        """ Export annotations from a model object to an SBML object

        Args:
            model_obj (:obj:`obj_model.Model`): model object
            nested_attr_paths (:obj:`dict`): dictionary that maps names of attributes to paths to the attributes
            sbml_obj (:obj:`libsbml.SBase`): SBML object

        Returns:
            :obj:`str`: string representation of the XML annotations of the model
        """
        key_vals = []
        for nested_attr_name, nested_attr_path in nested_attr_paths.items():
            attr = model_obj.get_nested_attr(nested_attr_path)
            val = model_obj.get_nested_attr_val(nested_attr_path)

            serialized_val = attr.serialize(val)

            if serialized_val is None:
                serialized_val = ''
            key_vals.append(('<wcLang:property>'
                             '<wcLang:key>{}</wcLang:key>'
                             '<wcLang:value>{}</wcLang:value>'
                             '</wcLang:property>').format(
                nested_attr_name, serialized_val))

        return ''.join(key_vals)

    @classmethod
    def get_annotations(cls, model_obj, nested_attr_paths, sbml_obj, model_objs=None):
        """ Import annotations from a model object to an SBML object

        Args:
            model_obj (:obj:`obj_model.Model`): model object
            nested_attr_paths (:obj:`dict`): dictionary that maps names of attributes to paths to the attributes
            sbml_obj (:obj:`libsbml.SBase`): SBML object
            model_objs (:obj:`dict`, optional): dictionary that maps classes of model objects to dictonaries
                that map ids of model objects to model objects
        """
        for nested_attr_path in nested_attr_paths.values():
            attr = model_obj.get_nested_attr(nested_attr_path)
            val = attr.get_none_value()
            model_obj.set_nested_attr_val(nested_attr_path, val)

        for nested_attr_name, val in cls.parse_annotations(sbml_obj).items():
            nested_attr_path = nested_attr_paths.get(nested_attr_name, None)
            if nested_attr_path:
                attr = model_obj.get_nested_attr(nested_attr_path)
                if isinstance(attr, obj_model.RelatedAttribute):
                    val, error = attr.deserialize(val, model_objs)
                else:
                    val, error = attr.deserialize(val)
                assert error is None, 'Error parsing {}.{} from {}: {}'.format(
                    model_obj.__class__.__name__, attr_name, str(val), str(error))
                model_obj.set_nested_attr_val(nested_attr_path, val)

    @classmethod
    def parse_annotations(cls, sbml_obj):
        """ Import annotations from a model object to an SBML object

        Args:
            sbml_obj (:obj:`libsbml.SBase`): SBML object

        Returns:
            :obj:`dict`: dictionary of names and values of annotated attributes
        """
        sbml_annots = cls.call_libsbml(sbml_obj.getAnnotation)
        key_vals = {}

        for i_annot in range(cls.call_libsbml(sbml_annots.getNumChildren, returns_int=True)):
            sbml_annot = cls.call_libsbml(sbml_annots.getChild, i_annot)
            if cls.call_libsbml(sbml_annot.getName) == 'annotation' and cls.call_libsbml(sbml_annot.getURI) == cls.XML_NAMESPACE:

                for i_child in range(cls.call_libsbml(sbml_annot.getNumChildren, returns_int=True)):
                    key_val = cls.call_libsbml(sbml_annot.getChild, i_child)

                    for i_g_child in range(cls.call_libsbml(key_val.getNumChildren, returns_int=True)):
                        g_child = cls.call_libsbml(key_val.getChild, i_g_child)
                        if cls.call_libsbml(g_child.getName) == 'key':
                            attr_name = cls.call_libsbml(g_child.toXMLString)[12:-13]
                        elif cls.call_libsbml(g_child.getName) == 'value':
                            val = cls.call_libsbml(g_child.toXMLString)[14:-15]

                    key_vals[attr_name] = val

        return key_vals

    @staticmethod
    def gen_nested_attr_paths(dotted_attr_paths):
        """ Generate structure representations of a list of nested attribute paths

        Args:
            dotted_attr_paths (:obj:`list` of :obj:`str`): List of paths to nested attributes. Each path is a
                dot-separated string of the series of nested attributes.

        Returns:
            :obj:`dict`: dictionary that ma
        """
        nested_attr_paths = {}
        for dotted_attr_path in dotted_attr_paths:
            nested_attr_paths[dotted_attr_path] = [(attr_name, ) for attr_name in dotted_attr_path.split('.')]
        return nested_attr_paths

    @classmethod
    def gen_authors_annotation(cls, model):
        """ Generate an annotation of the authors of a model

        Args:
            model (:obj:`wc_lang.core.Model`): model

        Returns:
            :obj:`str`: string representation of the XML annotation of the authors of a model
        """
        authors_sbml = ''
        for au in model.authors:
            authors_sbml += ('<wcLang:author>'
                             '<wcLang:id>{}</wcLang:id>'
                             '<wcLang:name>{}</wcLang:name>'
                             '<wcLang:lastName>{}</wcLang:lastName>'
                             '<wcLang:firstName>{}</wcLang:firstName>'
                             '<wcLang:middleName>{}</wcLang:middleName>'
                             '<wcLang:title>{}</wcLang:title>'
                             '<wcLang:organization>{}</wcLang:organization>'
                             '<wcLang:email>{}</wcLang:email>'
                             '<wcLang:website>{}</wcLang:website>'
                             '<wcLang:address>{}</wcLang:address>'
                             '<wcLang:identifiers>{}</wcLang:identifiers>'
                             '<wcLang:comments>{}</wcLang:comments>'
                             '</wcLang:author>').format(
                au.id, au.name, au.last_name, au.first_name, au.middle_name,
                au.title, au.organization, au.email, au.website, au.address,
                wc_lang.core.Author.Meta.attributes['identifiers'].serialize(au.identifiers),
                au.comments)
        return '<wcLang:authors>{}</wcLang:authors>'.format(authors_sbml)

    @classmethod
    def get_authors_annotation(cls, model, sbml_model, model_objs):
        """ Import the SBML annotation of the authors of a model to a model

        Args:
            model (:obj:`wc_lang.core.Model`): model
            sbml_model (:obj:`libsbml.Model`): SBML model
            model_objs (:obj:`dict`, optional): dictionary that maps classes of model objects to dictonaries
                that map ids of model objects to model objects
        """
        sbml_annots = cls.call_libsbml(sbml_model.getAnnotation)
        for i_annot in range(cls.call_libsbml(sbml_annots.getNumChildren, returns_int=True)):
            sbml_annot = cls.call_libsbml(sbml_annots.getChild, i_annot)
            if cls.call_libsbml(sbml_annot.getName) == 'annotation' and cls.call_libsbml(sbml_annot.getURI) == cls.XML_NAMESPACE:
                for i_annot_prop in range(cls.call_libsbml(sbml_annot.getNumChildren, returns_int=True)):
                    sbml_authors = cls.call_libsbml(sbml_annot.getChild, i_annot_prop)
                    if cls.call_libsbml(sbml_authors.getName) == 'authors':
                        for i_author in range(cls.call_libsbml(sbml_authors.getNumChildren, returns_int=True)):
                            sbml_au = cls.call_libsbml(sbml_authors.getChild, i_author)
                            au = model.authors.create()
                            for i_prop in range(cls.call_libsbml(sbml_au.getNumChildren, returns_int=True)):
                                sbml_prop = cls.call_libsbml(sbml_au.getChild, i_prop)

                                key = cls.call_libsbml(sbml_prop.getName)

                                val = cls.call_libsbml(sbml_prop.toXMLString)
                                val = val.partition('>')[2]
                                val = val[0:val.rfind('<')]

                                if key == 'id':
                                    au.id = val
                                elif key == 'name':
                                    au.name = val
                                elif key == 'lastName':
                                    au.last_name = val
                                elif key == 'firstName':
                                    au.first_name = val
                                elif key == 'middleName':
                                    au.middle_name = val
                                elif key == 'title':
                                    au.title = val
                                elif key == 'organization':
                                    au.organization = val
                                elif key == 'email':
                                    au.email = val
                                elif key == 'website':
                                    au.website = val
                                elif key == 'address':
                                    au.address = val
                                elif key == 'identifiers':
                                    attr = wc_lang.core.Author.Meta.attributes['identifiers']
                                    au.identifiers, error = attr.deserialize(val, model_objs)
                                    assert error is None, str(error)
                                elif key == 'comments':
                                    au.comments = val

    @classmethod
    def set_commments(cls, model_obj, sbml_obj):
        """ Export comments from a model object to an SBML object

        Args:
            model_obj (:obj:`obj_model.Model`): model object
            sbml_obj (:obj:`libsbml.SBase`): SBML object
        """
        if model_obj.comments:
            cls.call_libsbml(sbml_obj.setNotes, cls.str_to_xml_node(model_obj.comments))

    @classmethod
    def get_commments(cls, model_obj, sbml_obj):
        """ Import comments from an SBML object to a model object

        Args:
            model_obj (:obj:`obj_model.Model`): model object
            sbml_obj (:obj:`libsbml.SBase`): SBML object
        """
        if cls.call_libsbml(sbml_obj.isSetNotes):
            sbml_notes = cls.call_libsbml(sbml_obj.getNotes)
            model_obj.comments = cls.str_from_xml_node(cls.call_libsbml(sbml_notes.getChild, 0))
        else:
            model_obj.comments = ''

    @classmethod
    def call_libsbml(cls, method, *args, returns_int=False, debug=False):
        """ Call a libSBML method and handle any errors.

        Unfortunately, libSBML methods that do not return data usually report errors via return codes,
        instead of exceptions, and the generic return codes contain virtually no information.
        This function wraps these methods and raises useful exceptions when errors occur.

        Set `returns_int` `True` to avoid raising false exceptions or warnings from methods that return
        integer values.

        Args:
            method (:obj:`type`, :obj:`types.FunctionType`, or :obj:`types.MethodType`): `libsbml` method to execute
            args (:obj:`list`): a `list` of arguments to the `libsbml` method
            returns_int (:obj:`bool`, optional): whether the method returns an integer; if `returns_int`
                is `True`, then an exception will not be raised if the method call returns an integer
            debug (:obj:`bool`, optional): whether to print debug output

        Returns:
            :obj:`obj` or `int`: if the call does not return an error, return the `libsbml`
                method's return value, either an object that has been created or retrieved, or an integer
                value, or the `libsbml` success return code, :obj:`libsbml.LIBSBML_OPERATION_SUCCESS`

        Raises:
            :obj:`LibSbmlError`: if the `libsbml` call raises an exception, or returns None, or
                returns a known integer error code != :obj:`libsbml.LIBSBML_OPERATION_SUCCESS`
        """
        new_args = []
        for arg in args:
            new_args.append(arg)
        if new_args:
            new_args_str = ', '.join([str(a) for a in new_args])
            call_str = "method: {}; args: {}".format(method, new_args_str)
        else:
            call_str = "method: {}".format(method)
        if debug:
            print('libSBML call:', call_str)
        try:
            rc = method(*tuple(new_args))
        except BaseException as error:
            raise LibSbmlError("Error '{}' in libSBML method call '{}'.".format(error, call_str))
        if rc == None:
            raise LibSbmlError("libSBML returned None when executing '{}'.".format(call_str))
        elif type(rc) is int:
            # if `method` returns an int value, do not interpret rc as an error code
            if returns_int:
                if debug:
                    print('libSBML returns an int:', rc)
                return rc

            if rc == libsbml.LIBSBML_OPERATION_SUCCESS:
                if debug:
                    print('libSBML returns: LIBSBML_OPERATION_SUCCESS')
                return rc
            else:
                error_code = libsbml.OperationReturnValue_toString(rc)
                if error_code is None:
                    if debug:
                        print("libSBML returns:", rc)
                    warnings.warn("call_libsbml: unknown error code {} returned by '{}'."
                                  "\nPerhaps an integer value is being returned; if so, to avoid this warning "
                                  "pass 'returns_int=True' to call_libsbml().".format(error_code, call_str),
                                  wc_lang.core.WcLangWarning)
                    return rc
                else:
                    raise LibSbmlError("LibSBML returned error code '{}' when executing '{}'."
                                       "\nWARNING: if this libSBML call returns an int value, then this error may be "
                                       "incorrect; to avoid this error pass 'returns_int=True' to call_libsbml().".format(
                                           error_code, call_str))
        else:
            # return data provided by libSBML method
            if debug:
                print('libSBML returns:', rc)
            return rc

    @classmethod
    def str_to_xml_node(cls, str):
        """ Convert a Python string to an XML string that can be stored as a Note in an SBML document.

        Args:
            str (:obj:`str`): a string

        Returns:
            :obj:`libsbml.XMLNode`: an XML string that can be stored as a `Note` in an SBML document
        """
        return cls.call_libsbml(libsbml.XMLNode.convertStringToXMLNode,
                                '<p xmlns="http://www.w3.org/1999/xhtml">{}</p>'.format(str))

    @classmethod
    def str_from_xml_node(cls, xml_node):
        """ Convert an XML string (e.g., from a Note in an SBML document) to a Python string.

        Args:
            xml_node (:obj:`libsbml.XMLNode`): an XML string that can be stored as a `Note` in an SBML document

        Returns:
            :obj:`str`: a string
        """
        prefix = '<p xmlns="http://www.w3.org/1999/xhtml">'
        suffix = '</p>'
        text = cls.call_libsbml(xml_node.toXMLString)
        text = text[len(prefix):-len(suffix)]
        text = text.replace('\n  ', '\n').strip()
        return text

    @classmethod
    def raise_if_error(cls, sbml_doc, message):
        """ Raise an error, if an SBML object has errors

        Args:
            sbml_doc (:obj:`libsbml.SBMLDocument`): SBML document
            message (:obj:`str`): summary of error

        Raises:
            :obj:`LibSbmlError`: if the SBML object has errors
        """
        errors, warns, log_errors, log_warns = cls.get_errors_warnings(sbml_doc)

        if warns or log_warns:
            warnings.warn('{}:{}{}'.format(message,
                                           sbml_doc.__class__.__name__,
                                           ''.join(warns),
                                           ''.join(log_warns)),
                          wc_lang.core.WcLangWarning)

        if errors or log_errors:
            raise LibSbmlError('{}:{}{}'.format(message,
                                                sbml_doc.__class__.__name__,
                                                ''.join(errors),
                                                ''.join(log_errors)))

    @classmethod
    def get_errors_warnings(cls, sbml_doc):
        """ Get libSBML errors and warnings

        Args:
            sbml_doc (:obj:`libsbml.SBMLDocument`): SBML document

        Returns:
            :obj:`list` of :obj:`str`: error messages
            :obj:`list` of :obj:`str`: warning messages
            :obj:`list` of :obj:`str`: log error messages
            :obj:`list` of :obj:`str`: log warning messages
        """
        errors = []
        warns = []
        n_errors = cls.call_libsbml(sbml_doc.getNumErrors, returns_int=True)
        for i_error in range(n_errors):
            error = cls.call_libsbml(sbml_doc.getError, i_error)
            msg = '\n  {}: {}: {}'.format(cls.call_libsbml(error.getSeverityAsString),
                                          cls.call_libsbml(error.getShortMessage),
                                          cls.call_libsbml(error.getMessage))
            if cls.call_libsbml(error.getSeverity, returns_int=True) in [libsbml.LIBSBML_SEV_INFO, libsbml.LIBSBML_SEV_WARNING]:
                warns.append(msg)
            else:
                errors.append(msg)

        error_log = cls.call_libsbml(sbml_doc.getErrorLog)
        n_log_errors = cls.call_libsbml(error_log.getNumErrors, returns_int=True)
        log_errors = []
        log_warns = []
        for i_log_error in range(n_log_errors):
            log_error = cls.call_libsbml(error_log.getError, i_log_error)
            msg = '\n  {}: {}: {}'.format(cls.call_libsbml(log_error.getSeverityAsString),
                                          cls.call_libsbml(log_error.getShortMessage),
                                          cls.call_libsbml(log_error.getMessage))
            if cls.call_libsbml(log_error.getSeverity, returns_int=True) in [libsbml.LIBSBML_SEV_INFO, libsbml.LIBSBML_SEV_WARNING]:
                log_warns.append(msg)
            else:
                log_errors.append(msg)

        return (errors, warns, log_errors, log_warns)
