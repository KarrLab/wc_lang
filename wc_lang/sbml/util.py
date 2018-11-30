""" Utilities for writing a `wc_lang` model to SBML

Includes

* Exception definitions for `wc_lang.sbml`
* Higher level functions for creating SBML objects
* Utilities for wrapping libSBML calls and initializing libSBML models

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-11-01
:Copyright: 2017, Karr Lab
:License: MIT
"""

from libsbml import (LIBSBML_OPERATION_SUCCESS, OperationReturnValue_toString,
                     SBMLNamespaces, SBMLDocument)
from warnings import warn
import libsbml
import six
import wc_lang.core

# Centralize code that depends on levels and versions
# SBML level and version
SBML_LEVEL = 3
SBML_VERSION = 1
# Flux balance constraints (fbc) version
FBC_VERSION = 2

# SBML compatibility method for the version being used


def get_SBML_compatibility_method(sbml_document):
    return sbml_document.checkL3v1Compatibility


class Error(Exception):
    '''Base class for `libSBML` exceptions
    '''
    pass


class LibSBMLError(Error):
    '''Exception raised when libSBML returns an error
    '''

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        '''Provide the Exception's msg; needed for Python 2.7, although not documented
        '''
        return self.msg


class LibSBMLInterface(object):
    '''Methods that compactly use libSBML to create SBML objects.

    The libSBML method calls provide horribly narrow interfaces, typically exchanging one
    value per call, which creates extremely verbose code. These methods aggregate multiple
    libSBML method calls to enable more compact usage.
    '''

    @staticmethod
    def _create_sbml_doc_w_fbc():
        """ Create an SBMLDocument that uses the 'Flux Balance Constraints' extension.

        Returns:
            :obj:`libsbml.Unit`: the new SBML Document

        Raises:
            :obj:`LibSBMLError`: if a libSBML calls fails
        """
        sbmlns = wrap_libsbml(SBMLNamespaces, SBML_LEVEL, SBML_VERSION, 'fbc', FBC_VERSION)
        sbml_document = wrap_libsbml(SBMLDocument, sbmlns)
        wrap_libsbml(sbml_document.setPackageRequired, 'fbc', False)
        return sbml_document

    @staticmethod
    def _add_sbml_unit(unit_definition, unit_kind, exponent=1, scale=0, multiplier=1.0):
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
            :obj:`LibSBMLError`: if a libSBML calls fails
        """
        unit = wrap_libsbml(unit_definition.createUnit)
        wrap_libsbml(unit.setKind, unit_kind)
        wrap_libsbml(unit.setExponent, exponent)
        wrap_libsbml(unit.setScale, scale)
        wrap_libsbml(unit.setMultiplier, multiplier)
        return unit

    @staticmethod
    def _create_sbml_parameter(sbml_model, id, value, units, name=None, constant=True):
        """ Add an SBML Parameter to an SBML model.

        See the `libSBML documentation
        <http://sbml.org/Software/libSBML/docs/python-api/classlibsbml_1_1_parameter.html/>`_
        and the SBML specs for details.

        Args:
            sbml_model (:obj:`libsbml.Model`): a libSBML Model
            id (:obj:`str`): the id of the new SBML Parameter
            value (:obj:`obj`): the value of the new SBML Parameter
            units (:obj:`str`): the units of the new SBML Parameter
            name (:obj:`str`, optional): the name of the new SBML Parameter
            constant (:obj:`str`, optional): whether the new SBML Parameter is a constant

        Returns:
            :obj:`libsbml.Parameter`: the new SBML Parameter

        Raises:
            :obj:`LibSBMLError`: if a libSBML calls fails
            :obj:`ValueError`: if a `Parameter` with id `id` is already in use
        """
        try:
            wrap_libsbml(sbml_model.getParameter, id)
            raise ValueError("warning: '{}' is already in use as a Parameter id.".format(id))
        except LibSBMLError as e:
            sbml_parameter = wrap_libsbml(sbml_model.createParameter)
            wrap_libsbml(sbml_parameter.setIdAttribute, id)
            if not name is None:
                wrap_libsbml(sbml_parameter.setName, name)
            if not value is None:
                wrap_libsbml(sbml_parameter.setValue, value)
            if not units is None:
                wrap_libsbml(sbml_parameter.setUnits, units)
            wrap_libsbml(sbml_parameter.setConstant, constant)
            return sbml_parameter


create_sbml_doc_w_fbc = LibSBMLInterface._create_sbml_doc_w_fbc
add_sbml_unit = LibSBMLInterface._add_sbml_unit
create_sbml_parameter = LibSBMLInterface._create_sbml_parameter


def wrap_libsbml(method, *args, **kwargs):
    """ Wrap a libSBML method so that errors in return code can be easily handled.

    Unfortunately, libSBML methods that do not return data usually report errors via return codes,
    instead of exceptions, and the generic return codes contain virtually no information.
    This function wraps these methods and raises useful exceptions when errors occur.

    Set `returns_int` `True` to avoid raising false exceptions or warnings from methods that return
    integer values.

    Args:
        method (:obj:`obj`): a reference to the `libsbml` method to execute
        args (:obj:`list` of :obj:`obj`): a `list` of arguments to the `libsbml` method
        kwargs (:obj:`dict` of `obj`): a `dict` of options:
        returns_int (:obj:`bool`, optional): whether the method returns an integer; if `returns_int`
            is `True`, then an exception will not be raised if the method call returns an integer
        debug (:obj:`bool`, optional): whether to print debug output

    Returns:
        :obj:`obj` or `int`: if the call does not return an error, return the `libsbml`
        method's return value, either an object that has been created or retrieved, or an integer
        value, or the `libsbml` success return code, `LIBSBML_OPERATION_SUCCESS`

    Raises:
        :obj:`LibSBMLError`: if the `libsbml` call raises an exception, or returns None, or
        returns a known integer error code != `LIBSBML_OPERATION_SUCCESS`
    """
    # process kwargs
    returns_int = False
    if 'returns_int' in kwargs:
        returns_int = kwargs['returns_int']
        del kwargs['returns_int']
    debug = False
    if 'debug' in kwargs:
        debug = kwargs['debug']
        del kwargs['debug']

    # warn about unused kwargs
    for k in kwargs.keys():
        warn("wrap_libsbml: unknown kwargs key '{}'".format(k))

    new_args = []
    for arg in args:
        # if on Python 2, convert unicode text to str(), because libSBML doesn't use SWIG_PYTHON_2_UNICODE
        if six.PY2 and isinstance(arg, six.text_type):
            new_args.append(str(arg))  # pragma: no cover # Python 2 only
        else:
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
        raise LibSBMLError("Error '{}' in libSBML method call '{}'.".format(error, call_str))
    if rc == None:
        raise LibSBMLError("libSBML returned None when executing '{}'.".format(call_str))
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
                warn("wrap_libsbml: unknown error code {} returned by '{}'."
                     "\nPerhaps an integer value is being returned; if so, to avoid this warning "
                     "pass 'returns_int=True' to wrap_libsbml().".format(error_code, call_str))
                return rc
            else:
                raise LibSBMLError("LibSBML returned error code '{}' when executing '{}'."
                                   "\nWARNING: if this libSBML call returns an int value, then this error may be "
                                   "incorrect; to avoid this error pass 'returns_int=True' to wrap_libsbml().".format(
                                       error_code, call_str))
    else:
        # return data provided by libSBML method
        if debug:
            print('libSBML returns:', rc)
        return rc


def init_sbml_model(sbml_document):
    """ Create and initialize an SMBL model.

    Args:
         sbml_document (:obj:`obj`): a `libsbml` SBMLDocument

    Returns:
        :obj:`libsbml.model`: the SBML model

    Raises:
        :obj:`LibSBMLError`: if calling `libsbml` raises an error
    """
    # Modified copy of libsbml-5.15.0/examples/python/createSimpleModel.py from 2017-10-02
    sbml_model = wrap_libsbml(sbml_document.createModel)
    fbc_model_plugin = wrap_libsbml(sbml_model.getPlugin, 'fbc')
    wrap_libsbml(fbc_model_plugin.setStrict, True)

    # To produce a model with complete units for the reaction rates, we need
    # to set the 'timeUnits' and 'extentUnits' attributes on Model.  We
    # set 'substanceUnits' too, for good measure, though it's not strictly
    # necessary here because we also set the units for invididual species
    # in their definitions.
    wrap_libsbml(sbml_model.setTimeUnits, 'second')
    wrap_libsbml(sbml_model.setExtentUnits, 'mole')
    wrap_libsbml(sbml_model.setSubstanceUnits, 'mole')
    wrap_libsbml(sbml_model.setVolumeUnits, 'litre')

    # Create a unit definition we will need later.  Note that SBML Unit
    # objects must have all four attributes 'kind', 'exponent', 'scale'
    # and 'multiplier' defined.
    per_second = wrap_libsbml(sbml_model.createUnitDefinition)
    wrap_libsbml(per_second.setIdAttribute, 'per_second')
    add_sbml_unit(per_second, libsbml.UNIT_KIND_SECOND, exponent=-1)

    for unit_def, unit_def_meta in wc_lang.core.ConcentrationUnit.Meta.items():
        sbml_unit_def = wrap_libsbml(sbml_model.createUnitDefinition)
        wrap_libsbml(sbml_unit_def.setIdAttribute, unit_def_meta['xml_id'])
        substance_unit = unit_def_meta['substance_units']
        add_sbml_unit(sbml_unit_def,
                      getattr(libsbml, 'UNIT_KIND_' + substance_unit['kind'].upper()),
                      exponent=substance_unit['exponent'],
                      scale=substance_unit['scale'])

    mmol_per_gDW_per_hr = wrap_libsbml(sbml_model.createUnitDefinition)
    wrap_libsbml(mmol_per_gDW_per_hr.setIdAttribute, 'mmol_per_gDW_per_hr')
    add_sbml_unit(mmol_per_gDW_per_hr, libsbml.UNIT_KIND_MOLE, scale=-3)
    add_sbml_unit(mmol_per_gDW_per_hr, libsbml.UNIT_KIND_GRAM, exponent=-1)
    add_sbml_unit(mmol_per_gDW_per_hr, libsbml.UNIT_KIND_SECOND, exponent=-1,
                  multiplier=3600.0)

    dimensionless = wrap_libsbml(sbml_model.createUnitDefinition)
    wrap_libsbml(dimensionless.setIdAttribute, 'dimensionless_ud')
    add_sbml_unit(dimensionless, libsbml.UNIT_KIND_DIMENSIONLESS)

    return sbml_model


def str_to_xmlstr(str):
    """ Convert a Python string to an XML string that can be stored as a Note in an SBML Document.

    Args:
        str (:obj:`str`): a string

    Returns:
        :obj:`str`: an XML string that can be stored as a Note in an SBML Document
    """
    # TODO: GET libSBML to do this XML crap, but none of the obvious methods work
    return "<p xmlns=\"http://www.w3.org/1999/xhtml\">{}</p>".format(str)
