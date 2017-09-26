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


def __wrap_libsbml(call, globals, locals):
    """ Wrap a libsbml method and properly handle errors.

    Unfortunately, libsbml methods that do not return data usually handle errors via return codes,
    instead of exceptions, and the generic return codes contain virtually no information.
    This function wraps these methods and raises useful exceptions when errors occur.

    Args:
        call (:obj:`str`): a libsbml expression to execute
        globals (:obj:`namespace`): the global namespace
        locals (:obj:`namespace`): the local namespace at the calling code

    Returns:
        :obj:`obj` or `int`: return the value returned by the libsbml method, either
        an object that has been created or retrieved, or an integer return code

    Raises:
        :obj:`LibSBMLError`: if `call` contains an error, or the libsbml call returns None,
        or the libsbml call return a code != LIBSBML_OPERATION_SUCCESS
    """
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
            return rc
        else:
            raise LibSBMLError("LibSBML returned error code '{}' "
                "when executing '{}'.".format(OperationReturnValue_toString(rc), call))
    else:
        # return data provided by libsbml method
        return rc

def wrap_libsbml(call):
    """ Wrap a libsbml method, automatically passing global and local namespaces.

    Args:
        call (:obj:`str`): the libsbml expression to execute

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
            frame.f_back.f_locals)
    finally:
        del frame

def str_to_xmlstr(str):
    """ Convert a Python string to an XML string that can be stored as a Note in an SBML Document.

    Args:
        str (:obj:`str`): a string

    Returns:
        :obj:`str`: an XML string that can be stored as a Note in an SBML Document
    """
    # TODO: GET libsbml to do this XML crap, but none of the obvious methods work
    return "<p xmlns=\"http://www.w3.org/1999/xhtml\">{}</p>".format(str)
