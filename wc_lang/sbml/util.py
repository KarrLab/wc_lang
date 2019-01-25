""" Utilities for writing/reading a `wc_lang` model to/from SBML

Includes

* Exception definitions for `wc_lang.sbml`
* Higher level functions for creating SBML objects
* Utilities for wrapping libSBML calls and initializing libSBML models

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-24
:Copyright: 2017-2019, Karr Lab
:License: MIT
"""

from libsbml import (LIBSBML_OPERATION_SUCCESS, OperationReturnValue_toString,
                     SBMLNamespaces, SBMLDocument)
import libsbml
import math
import pint
import types
import warnings
import wc_lang.core


class LibSbmlError(Exception):
    ''' Exception raised when libSBML returns an error '''

    def __init__(self, msg):
        self.msg = msg


class LibSbmlInterface(object):
    '''Methods for compactly using libSBML to create SBML objects.

    The libSBML method calls provide narrow interfaces, typically exchanging one
    value per call, which creates verbose code. The methods below aggregate multiple
    libSBML method calls to enable more compact usage.
    '''

    @classmethod
    def create_doc(cls, level=3, version=1, packages=None):
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
        namespaces = [SBMLNamespaces, level, version]
        for package_id, package_version in packages.items():
            namespaces.append(package_id)
            namespaces.append(package_version)
        sbml_ns = cls.call_libsbml(*namespaces)

        # create SBML document
        doc = cls.call_libsbml(SBMLDocument, sbml_ns)

        # set package requirements
        for package in packages:
            cls.call_libsbml(doc.setPackageRequired, package, False)

        # return SBML document
        return doc

    @classmethod
    def is_doc_compatible(cls, doc, level=3, version=1):
        """ Check the compatibility of an SBML document with a specific level and version

        Args:
            doc (:obj:`libsbml.SBMLDocument`): SBML document
            level (:obj:`int`, optional): SBML level number
            version (:obj:`int`, optional): SBML version number
        """
        # SBML compatibility method for the version being used
        return 0 == cls.call_libsbml(doc.checkL3v1Compatibility, returns_int=True)

    @classmethod
    def create_model(cls, doc):
        """ Create a SBML model

        Args:
            doc (:obj:`libsbml.SBMLDocument`): SBML document

        Returns:
            :obj:`libsbml.Model`: SBML model
        """
        return cls.call_libsbml(doc.createModel)

    @classmethod
    def init_model(cls, model, doc, packages=None):
        """ Create and initialize an SMBL model.

        Args:
            model (:obj:`wc_lang.core.Model`): model
            doc (:obj:`libsbml.SBMLDocument`): a `libsbml` SBMLDocument
            packages (:obj:`dict` that maps :obj:`str` to :obj:`int`, optional): dictionary of required packages
                that maps package identifiers to package numbers

        Returns:
            :obj:`libsbml.Model`: the SBML model

        Raises:
            :obj:`ValueError`: if unit is not supported
            :obj:`LibSbmlError`: if calling `libsbml` raises an error
        """

        # create model
        sbml_model = cls.create_model(doc)

        # enable plugins for packages
        packages = packages or {}
        for package_id in packages.keys():
            plugin = cls.call_libsbml(sbml_model.getPlugin, package_id)
            cls.call_libsbml(plugin.setStrict, True)

        # Set units
        cls.set_units(model, sbml_model)

        # return SBML model
        return sbml_model

    @classmethod
    def set_units(cls, model, sbml_model):
        """ Set time, extent, and substance units

        Args:
            model (:obj:`wc_lang.core.Model`): model
            sbml_model (:obj:`libsbml.Model`): SBML model that encodes the model
        """
        cls.call_libsbml(sbml_model.setTimeUnits, 'second')
        cls.call_libsbml(sbml_model.setExtentUnits, 'mole')
        cls.call_libsbml(sbml_model.setSubstanceUnits, 'mole')
        cls.call_libsbml(sbml_model.setVolumeUnits, 'litre')

        # Define units
        # Note that SBML Unit objects must have four attributes:
        # * `kind`
        # * `exponent`
        # * `scale`
        # * `multiplier`
        per_second = cls.call_libsbml(sbml_model.createUnitDefinition)
        cls.call_libsbml(per_second.setIdAttribute, 'per_second')
        cls.add_unit(per_second, libsbml.UNIT_KIND_SECOND, exponent=-1)

        unit_registry = wc_lang.core.DistributionInitConcentration.units.registry
        molar = unit_registry.parse_expression('M')
        molecule = unit_registry.parse_expression('molecule')
        for unit in wc_lang.core.DistributionInitConcentration.units.choices:
            sbml_unit_def = cls.call_libsbml(sbml_model.createUnitDefinition)
            cls.call_libsbml(sbml_unit_def.setIdAttribute, str(unit))

            if not isinstance(unit, unit_registry.Unit):
                raise ValueError('Unsupported units "{}"'.format(
                    unit))  # pragma: no cover # unreachable because all choices in above two cases
            expr = unit_registry.parse_expression(str(unit))

            try:
                molar_units = unit_registry.parse_units('molar')
                conversion = expr.to(molar_units)
                scale = math.log10(conversion.magnitude)
                assert int(scale) == scale
                cls.add_unit(sbml_unit_def,
                             getattr(libsbml, 'UNIT_KIND_MOLE'),
                             exponent=1.,
                             scale=int(scale))
                continue
            except pint.DimensionalityError:
                pass

            try:
                molecule_units = unit_registry.parse_units('molecule')
                conversion = expr.to(molecule_units)
                scale = math.log10(conversion.magnitude)
                assert int(scale) == scale
                cls.add_unit(sbml_unit_def,
                             getattr(libsbml, 'UNIT_KIND_ITEM'),
                             exponent=1.,
                             scale=int(scale))
                continue
            except pint.DimensionalityError:
                pass

            raise ValueError('Invalid unit {}'.format(str(unit)))  # pragma: no cover # unreachable because all choices in above two cases

        mmol_per_gDW_per_hr = cls.call_libsbml(sbml_model.createUnitDefinition)
        cls.call_libsbml(mmol_per_gDW_per_hr.setIdAttribute, 'mmol_per_gDW_per_hr')
        cls.add_unit(mmol_per_gDW_per_hr, libsbml.UNIT_KIND_MOLE, scale=-3)
        cls.add_unit(mmol_per_gDW_per_hr, libsbml.UNIT_KIND_GRAM, exponent=-1)
        cls.add_unit(mmol_per_gDW_per_hr, libsbml.UNIT_KIND_SECOND, exponent=-1,
                     multiplier=3600.0)

        dimensionless = cls.call_libsbml(sbml_model.createUnitDefinition)
        cls.call_libsbml(dimensionless.setIdAttribute, 'dimensionless_ud')
        cls.add_unit(dimensionless, libsbml.UNIT_KIND_DIMENSIONLESS)

    @classmethod
    def add_unit(cls, unit_definition, unit_kind, exponent=1, scale=0, multiplier=1.0):
        """ Add an SBML `Unit` to an existing SBML `UnitDefinition`.

        Provides the SBML level 3 version 1 default values for `exponent=1`, `scale=0`, and `multiplier=1.0`.
        See the `libSBML documentation
        <http://sbml.org/Software/libSBML/docs/python-api/classlibsbml_1_1_unit_definition.html/>`_
        and the SBML specs for details.

        Args:
            unit_definition (:obj:`libsbml.UnitDefinition`): a libSBML `UnitDefinition`
            unit_kind (:obj:`libsbml.UnitDefinition`): the unit kind code for the SBML unit
            exponent (:obj:`int`, optional): the exponent on the SBML unit
            scale (:obj:`int`, optional): the scale of the SBML unit
            multiplier (:obj:`float`, optional): the multiplier of the SBML unit

        Returns:
            :obj:`libsbml.Unit`: the new SBML Unit

        Raises:
            :obj:`LibSbmlError`: if a libSBML calls fails
        """
        unit = cls.call_libsbml(unit_definition.createUnit)
        cls.call_libsbml(unit.setKind, unit_kind)
        cls.call_libsbml(unit.setExponent, exponent)
        cls.call_libsbml(unit.setScale, scale)
        cls.call_libsbml(unit.setMultiplier, multiplier)
        return unit

    @classmethod
    def create_parameter(cls, model, id, value, units, name=None, constant=True):
        """ Add an SBML Parameter to an SBML model.

        See the `libSBML documentation
        <http://sbml.org/Software/libSBML/docs/python-api/classlibsbml_1_1_parameter.html/>`_
        and the SBML specs for details.

        Args:
            model (:obj:`libsbml.Model`): a libSBML Model
            id (:obj:`str`): the id of the new SBML Parameter
            value (:obj:`obj`): the value of the new SBML Parameter
            units (:obj:`str`): the units of the new SBML Parameter
            name (:obj:`str`, optional): the name of the new SBML Parameter
            constant (:obj:`str`, optional): whether the new SBML Parameter is a constant

        Returns:
            :obj:`libsbml.Parameter`: the new SBML Parameter

        Raises:
            :obj:`LibSbmlError`: if a libSBML calls fails
            :obj:`ValueError`: if a `Parameter` with id `id` is already in use
        """
        try:
            cls.call_libsbml(model.getParameter, id)
            raise ValueError("warning: '{}' is already in use as a Parameter id.".format(id))
        except LibSbmlError as e:
            sbml_parameter = cls.call_libsbml(model.createParameter)
            cls.call_libsbml(sbml_parameter.setIdAttribute, id)
            if not name is None:
                cls.call_libsbml(sbml_parameter.setName, name)
            if not value is None:
                cls.call_libsbml(sbml_parameter.setValue, value)
            if not units is None:
                cls.call_libsbml(sbml_parameter.setUnits, units)
            cls.call_libsbml(sbml_parameter.setConstant, constant)
            return sbml_parameter

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
                value, or the `libsbml` success return code, `LIBSBML_OPERATION_SUCCESS`

        Raises:
            :obj:`LibSbmlError`: if the `libsbml` call raises an exception, or returns None, or
            returns a known integer error code != `LIBSBML_OPERATION_SUCCESS`
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

            if rc == LIBSBML_OPERATION_SUCCESS:
                if debug:
                    print('libSBML returns: LIBSBML_OPERATION_SUCCESS')
                return rc
            else:
                error_code = OperationReturnValue_toString(rc)
                if error_code is None:
                    if debug:
                        print("libSBML returns:", rc)
                    warnings.warn("call_libsbml: unknown error code {} returned by '{}'."
                                  "\nPerhaps an integer value is being returned; if so, to avoid this warning "
                                  "pass 'returns_int=True' to call_libsbml().".format(error_code, call_str), UserWarning)
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

    @staticmethod
    def str_to_xmlstr(str):
        """ Convert a Python string to an XML string that can be stored as a Note in an SBML document.

        Args:
            str (:obj:`str`): a string

        Returns:
            :obj:`str`: an XML string that can be stored as a `Note` in an SBML document
        """
        # TODO: GET libSBML to do this XML crap, but none of the obvious methods work
        return "<p xmlns=\"http://www.w3.org/1999/xhtml\">{}</p>".format(str)
