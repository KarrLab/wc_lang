""" Utilities for writing/reading a wc_lang model to/from SBML

Representations include
* Files
* Strings

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

"""
License for code reused from libSBML:
<!--------------------------------------------------------------------------
This sample program is distributed under a different license than the rest
of libSBML.  This program uses the open-source MIT license, as follows:
##
Copyright (c) 2013-2017 by the California Institute of Technology
(California, USA), the European Bioinformatics Institute (EMBL-EBI, UK)
and the University of Heidelberg (Germany), with support from the National
Institutes of Health (USA) under grant R01GM070923.  All rights reserved.
##
Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:
##
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
##
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
##
Neither the name of the California Institute of Technology (Caltech), nor
of the European Bioinformatics Institute (EMBL-EBI), nor of the University
of Heidelberg, nor the names of any contributors, may be used to endorse
or promote products derived from this software without specific prior
written permission.
------------------------------------------------------------------------ -->
"""

import sys
import inspect
from libsbml import (LIBSBML_OPERATION_SUCCESS, UNIT_KIND_SECOND, UNIT_KIND_MOLE, UNIT_KIND_GRAM,
    UNIT_KIND_DIMENSIONLESS, OperationReturnValue_toString, SBMLNamespaces, SBMLDocument)

import six

# SBML level and version being used
SBML_LEVEL = 3
SBML_VERSION = 1

def get_SBML_compatibility_method(sbml_document):
    return sbml_document.checkL3v1Compatibility
    

class Error(Exception):
    '''Base class libsbml exceptions.'''
    pass


class LibSBMLError(Error):
    '''Exception raised when libsbml returns an error.'''

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        '''Provide the Exception's msg; needed for Python 2.7, although not documented.'''
        return self.msg

class LibsbmlInterface(object):
    '''Methods that compactly use libsbml to create SBML objects.

    The libsbml method calls provide horribly narrow interfaces, typically exchanging one
    value per call, which creates extremely verbose SBML code. These methods aggregate multiple
    libsbml method calls to enable more compact usage.
    '''

    @staticmethod
    def _create_sbml_doc_w_fbc():
        """ Create an SBMLDocument that uses the 'Flux Balance Constraints' extension, version 2.

        Returns:
            :obj:`libsbml.Unit`: the new SBML Document

        Raises:
            :obj:`LibSBMLError`: if one of the libsbml calls fails
        """
        sbmlns = wrap_libsbml(SBMLNamespaces, SBML_LEVEL, SBML_VERSION, 'fbc', 2)
        sbml_document = wrap_libsbml(SBMLDocument, sbmlns)
        wrap_libsbml(sbml_document.setPackageRequired, 'fbc', False)
        return sbml_document

    @staticmethod
    def _add_sbml_unit(unit_definition, unit_kind, exponent=1, scale=0, multiplier=1.0):
        """ Add an SBML unit on an existing SBML unit definition.

        Provides the SBML level 3 version 1 default values for `exponent=1`, `scale=0`, and `multiplier=1.0`.
        See http://sbml.org/Software/libSBML/docs/python-api/classlibsbml_1_1_unit_definition.html
        in the libsbml documentation and the SBML specs for details.

        Args:
            unit_definition (:obj:`libsbml.UnitDefinition`): a libsbml UnitDefinition
            unit_kind (:obj:`libsbml.UnitDefinition`): the unit kind code for the SBML unit
            exponent (:obj:`int`, optional): the exponent on the SBML unit
            scale (:obj:`int`, optional): the scale of the SBML unit
            multiplier (:obj:`float`, optional): the multiplier of the SBML unit

        Returns:
            :obj:`libsbml.Unit`: the new SBML Unit

        Raises:
            :obj:`LibSBMLError`: if one of the libsbml calls fails
        """
        unit = wrap_libsbml(unit_definition.createUnit)
        wrap_libsbml(unit.setKind, unit_kind)
        wrap_libsbml(unit.setExponent, exponent)
        wrap_libsbml(unit.setScale, scale)
        wrap_libsbml(unit.setMultiplier, multiplier)
        return unit

    @staticmethod
    def _create_sbml_parameter(sbml_model, id, name=None, value=None, units=None, constant=True):
        """ Add an SBML Parameter to an SBML model.

        See http://sbml.org/Software/libSBML/docs/python-api/classlibsbml_1_1_parameter.html
        in the libsbml documentation and the SBML specs for details.

        Args:
            sbml_model (:obj:`libsbml.Model`): a libsbml Model
            id (:obj:`str`): the id of the new SBML Parameter
            name (:obj:`str`, optional): the name of the new SBML Parameter
            value (:obj:`obj`, optional): the value of the new SBML Parameter
            units (:obj:`str`, optional): the units of the new SBML Parameter
            constant (:obj:`str`, optional): whether the new SBML Parameter is a constant

        Returns:
            :obj:`libsbml.Parameter`: the new SBML Parameter

        Raises:
            :obj:`LibSBMLError`: if one of the libsbml calls fails
            :obj:`ValueError`: if the Parameter `id` is already in use
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

create_sbml_doc_w_fbc = LibsbmlInterface._create_sbml_doc_w_fbc
add_sbml_unit = LibsbmlInterface._add_sbml_unit
create_sbml_parameter = LibsbmlInterface._create_sbml_parameter

def wrap_libsbml(method, *args, **kwargs):
    """ Wrap a libsbml method so that errors in return code can be easily handled.

    Unfortunately, libsbml methods that do not return data usually report errors via return codes,
    instead of exceptions, and the generic return codes contain virtually no information.
    This function wraps these methods and raises useful exceptions when errors occur.

    Args:
        method (:obj:`obj`): a reference to the `libsbml` method to execute
        args (:obj:`list` of `obj`): a `list` of arguments to the `libsbml` method
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
    debug = False
    if 'debug' in kwargs:
        debug = kwargs['debug']
    # TODO: report error on unused kwargs

    new_args = []
    for arg in args:
        # if on Python 2, convert unicode text to str(), because libsbml doesn't use SWIG right
        if six.PY2 and isinstance(arg, six.text_type):
            new_args.append(str(arg))
        else:
            new_args.append(arg)
    if debug:
        # TODO: make a string for arguments
        print('libsbml call:', method)
    try:
        rc = method(*tuple(args))
    except BaseException as error:
        raise LibSBMLError("Error '{}' in libsbml method call '{}'.".format(error, method))
    if rc == None:
        raise LibSBMLError("libsbml returned None when executing '{}'.".format(method))
    elif type(rc) is int:
        if rc == LIBSBML_OPERATION_SUCCESS:
            if debug:
                print('libsbml returns: LIBSBML_OPERATION_SUCCESS')
            return rc
        else:
            error_code = OperationReturnValue_toString(rc)
            # Handle libsbml methods that return int as values
            # TODO: handle this more gracefully
            if error_code is None or returns_int:
                if debug:
                    print("libsbml returns:", rc)
                return rc
            else:
                raise LibSBMLError("LibSBML returned error code '{}' "
                    "when executing '{}'.\nWARNING: if the libsbml call above returns an int, then this "
                    "error may be incorrect; pass 'returns_int=True' to wrap_libsbml().".format(error_code, method))
    else:
        # return data provided by libsbml method
        if debug:
            print('libsbml returns:', rc)
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
    # TODO: set strict to True
    wrap_libsbml(fbc_model_plugin.setStrict, False)

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
    add_sbml_unit(per_second, UNIT_KIND_SECOND, exponent=-1)

    mmol_per_gDW_per_hr = wrap_libsbml(sbml_model.createUnitDefinition)
    wrap_libsbml(mmol_per_gDW_per_hr.setIdAttribute, 'mmol_per_gDW_per_hr')
    add_sbml_unit(mmol_per_gDW_per_hr, UNIT_KIND_MOLE, scale=-3)
    add_sbml_unit(mmol_per_gDW_per_hr, UNIT_KIND_GRAM, exponent=-1)
    add_sbml_unit(mmol_per_gDW_per_hr, UNIT_KIND_SECOND, exponent=-1,
        multiplier=3600.0)

    dimensionless = wrap_libsbml(sbml_model.createUnitDefinition)
    wrap_libsbml(dimensionless.setIdAttribute, 'dimensionless_ud')
    add_sbml_unit(dimensionless, UNIT_KIND_DIMENSIONLESS)

    return sbml_model


def str_to_xmlstr(str):
    """ Convert a Python string to an XML string that can be stored as a Note in an SBML Document.

    Args:
        str (:obj:`str`): a string

    Returns:
        :obj:`str`: an XML string that can be stored as a Note in an SBML Document
    """
    # TODO: GET libsbml to do this XML crap, but none of the obvious methods work
    return "<p xmlns=\"http://www.w3.org/1999/xhtml\">{}</p>".format(str)
