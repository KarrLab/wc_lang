""" Reading/writing a WC model to/from an SBML representation

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
from libsbml import *

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
    def _create_sbml_unit(unit_definition, unit_kind, exponent=1, scale=0, multiplier=1.0):
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
        unit = wrap_libsbml("unit_definition.createUnit()")
        wrap_libsbml("unit.setKind({})".format(unit_kind))
        wrap_libsbml("unit.setExponent({})".format(exponent))
        wrap_libsbml("unit.setScale({})".format(scale))
        wrap_libsbml("unit.setMultiplier({})".format(multiplier))
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
            wrap_libsbml("sbml_model.getParameter('{}')".format(id))
            raise ValueError("warning: '{}' is already in use as a Parameter id.".format(id))
        except LibSBMLError as e:
            sbml_parameter = wrap_libsbml("sbml_model.createParameter()")
            wrap_libsbml("sbml_parameter.setIdAttribute('{}')".format(id))
            if not name is None:
                wrap_libsbml("sbml_parameter.setName('{}')".format(name))
            if not value is None:
                wrap_libsbml("sbml_parameter.setValue({})".format(value))
            if not units is None:
                wrap_libsbml("sbml_parameter.setUnits('{}')".format(units))
            wrap_libsbml("sbml_parameter.setConstant({})".format(constant))
            return sbml_parameter

create_sbml_unit = LibsbmlInterface._create_sbml_unit
create_sbml_parameter = LibsbmlInterface._create_sbml_parameter

def __wrap_libsbml(call, globals, locals, returns_int=False, debug=False):
    """ Wrap a libsbml method and properly handle errors.

    Unfortunately, libsbml methods that do not return data usually handle errors via return codes,
    instead of exceptions, and the generic return codes contain virtually no information.
    This function wraps these methods and raises useful exceptions when errors occur.

    Args:
        call (:obj:`str`): a libsbml expression to execute
        globals (:obj:`namespace`): the global namespace
        locals (:obj:`namespace`): the local namespace at the calling code
        returns_int (:obj:`bool`, optional): whether the method returns an int
        debug (:obj:`bool`, optional): whether to print debug output

    Returns:
        :obj:`obj` or `int`: return the value returned by the libsbml method, either
        an object that has been created or retrieved, or an integer return code

    Raises:
        :obj:`LibSBMLError`: if `call` contains an error, or the libsbml call returns None,
        or the libsbml call return a code != LIBSBML_OPERATION_SUCCESS
    """
    if debug:
        print('_wrap_libsbml call:', call)
    try:
        rc = eval(call, globals, locals)
    except SyntaxError as error:
        raise LibSBMLError("Syntax error in libsbml method call '{}'.".format(call))
    except NameError as error:
        raise LibSBMLError("NameError '{}' in libsbml method call '{}'.".format(error, call))
    except Exception as error:
        raise LibSBMLError("Error '{}' in libsbml method call '{}'.".format(error, call))
    if rc == None:
        raise LibSBMLError("libsbml returned None when executing '{}'.".format(call))
    elif type(rc) is int:
        if rc == LIBSBML_OPERATION_SUCCESS:
            if debug:
                print('libsbml returned LIBSBML_OPERATION_SUCCESS')
            return rc
        else:
            error_code = OperationReturnValue_toString(rc)
            # Handle libsbml methods that return int as values
            # TODO: handle this more gracefully
            if error_code is None or returns_int:
                if debug:
                    print("libsbml returns int that isn't error code", rc)
                return rc
            else:
                raise LibSBMLError("LibSBML returned error code '{}' "
                    "when executing '{}'.\nWARNING: if the libsbml call above returns an int, then this "
                    "error may be incorrect; pass 'returns_int=True' to wrap_libsbml().".format(error_code, call))
    else:
        # return data provided by libsbml method
        if debug:
            print('libsbml returned', rc)
        return rc

def wrap_libsbml(call, returns_int=False, debug=False):
    """ Wrap a libsbml method, automatically passing global and local namespaces.

    Args:
        call (:obj:`str`): the libsbml expression to execute
        returns_int (:obj:`bool`, optional): whether the method returns an int
        debug (:obj:`bool`, optional): whether to print debug output

    Returns:
        :obj:`obj` or `int`: return the libsbml method's return value, either
        an object that has been created or retrieved, or an integer return code

    Raises:
        :obj:`LibSBMLError`: if `call` contains an error, or the libsbml call returns None,
        or the libsbml call return a code != LIBSBML_OPERATION_SUCCESS
    """
    frame = inspect.currentframe()
    try:
        return __wrap_libsbml(call,
            frame.f_back.f_globals,
            frame.f_back.f_locals,
            returns_int, debug)
    finally:
        del frame

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
    sbml_model = wrap_libsbml("sbml_document.createModel()")

    # To produce a model with complete units for the reaction rates, we need
    # to set the 'timeUnits' and 'extentUnits' attributes on Model.  We
    # set 'substanceUnits' too, for good measure, though it's not strictly
    # necessary here because we also set the units for invididual species
    # in their definitions.
    wrap_libsbml("sbml_model.setTimeUnits('second')")
    wrap_libsbml("sbml_model.setExtentUnits('mole')")
    wrap_libsbml("sbml_model.setSubstanceUnits('mole')")
    wrap_libsbml("sbml_model.setVolumeUnits('litre')")

    # Create a unit definition we will need later.  Note that SBML Unit
    # objects must have all four attributes 'kind', 'exponent', 'scale'
    # and 'multiplier' defined.
    per_second = wrap_libsbml("sbml_model.createUnitDefinition()")
    wrap_libsbml("per_second.setId('per_second')")
    create_sbml_unit(per_second, UNIT_KIND_SECOND, exponent=-1)

    mmol_per_gDW_per_hr = wrap_libsbml("sbml_model.createUnitDefinition()")
    wrap_libsbml("mmol_per_gDW_per_hr.setId('mmol_per_gDW_per_hr')")
    create_sbml_unit(mmol_per_gDW_per_hr, UNIT_KIND_MOLE, scale=-3)
    create_sbml_unit(mmol_per_gDW_per_hr, UNIT_KIND_GRAM, exponent=-1)
    create_sbml_unit(mmol_per_gDW_per_hr, UNIT_KIND_SECOND, exponent=-1,
        multiplier=3600.0)

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
