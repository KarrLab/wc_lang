""" Reading/writing a WC model to/from an SBML representation

Representations include
* Files
* Strings

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import sys
from libsbml import *

from wc_lang.core import (Model, Taxon, Submodel, ObjectiveFunction, Compartment, SpeciesType,
    Species, Concentration, Reaction, ReactionParticipant, RateLaw, RateLawEquation,
    BiomassComponent, BiomassReaction, Parameter, Reference, CrossReference)

from wc_lang.sbml.util import wrap_libsbml

'''
    reader = SBMLReader()
    document = reader.readSBML(file)
    document.getNumErrors()
    writeSBMLToFile(SBMLDocument d, string filename)
'''

'''
wc_lang to SBML mapping to support FBA modeling
Individual wc_lang submodels that use dFBA are mapped to individual SBML documents and files.

WC			            SBML                                                Status
-----                   -----                                               ------
Model                   Model                                               Ignored
Taxon			        None, perhaps make SBML annotations                 Ignored
Submodel			    Model                                               Implemented
ObjectiveFunction       Objective                                           TBD
Compartment			    Compartment                                         Implemented
SpeciesType			    NA: SpeciesType aren't defined
Species			        Species                                             Implemented
Concentration			NA: concentrations are incorporated in Species
Reaction			    Reaction, but FbcReactionPlugin for DFBA submodels  TBD
ReactionParticipant		SpeciesReference in a Reaction
RateLaw			        KineticLaw?                                         Ignored
RateLawEquation			
BiomassComponent			
BiomassReaction			
Parameter			    Parameter
Reference			    
CrossReference			
'''

'''
wc_lang attribute to SBML mapping
WC Model			    SBML Model
comments                notes
references              notes, as a Python dict

'''

class Reader(object):
    """ Read model objects from an SBML representation """

    def run(self, objects, models, path=None, get_related=True,
        title=None, description=None, keywords=None, version=None,
        language=None, creator=None):
        """ Write a list of model objects to a string or a file in an SBML format.

        Args:
            objects (:obj:`list`): list of objects
            models (:obj:`list`): list of model, in the order that they should
                appear as worksheets; all models which are not in `models` will
                follow in alphabetical order
            path (:obj:`str`, optional): path of SBML file to write
            title (:obj:`str`, optional): title
            description (:obj:`str`, optional): description
            keywords (:obj:`str`, optional): keywords
            version (:obj:`str`, optional): version
            language (:obj:`str`, optional): language
            creator (:obj:`str`, optional): creator
        """
        pass

class Writer(object):
    """ Write model objects to an SBML representation """

    def run(self, objects, models, path=None, get_related=True,
        title=None, description=None, keywords=None, version=None,
        language=None, creator=None):
        """ Write a list of model objects to a string or a file in an SBML format.

        Args:
            objects (:obj:`list`): list of objects
            models (:obj:`list`): list of model, in the order that they should
                appear as worksheets; all models which are not in `models` will
                follow in alphabetical order
            path (:obj:`str`, optional): path of SBML file to write
            title (:obj:`str`, optional): title
            description (:obj:`str`, optional): description
            keywords (:obj:`str`, optional): keywords
            version (:obj:`str`, optional): version
            language (:obj:`str`, optional): language
            creator (:obj:`str`, optional): creator
        """
        pass

class SBMLExchange(object):
    """ Exchange `wc_lang` model to/from a libSBML SBML representation """

    @staticmethod
    def write(objects, models):
        """ Write a list of `wc_lang` model objects to a libSBML `SBMLDocument`.

        Warning: `wc_lang` and SBML semantics are not equivalent. Thus, this
        `wc_lang` information will not be written:
            * xxx

        Args:
            objects (:obj:`list`): list of objects
            models (:obj:`list`): list of `Model`

        Returns:
            :obj:`SBMLDocument`: an SBMLDocument containing `objects`

        Raises:
            :obj:`ValueError`: if ...
        """
        # initialize SBMLDocument
        # Create an empty SBMLDocument object.  It's a good idea to check for
        # possible errors.  Even when the parameter values are hardwired like
        # this, it is still possible for a failure to occur (e.g., if the
        # operating system runs out of memory).
        try:
            document = SBMLDocument(3, 2)
        except ValueError:
            raise SystemExit('Could not create SBMLDocumention object')

        # Create the basic Model object inside the SBMLDocument object.
        sbml_model = wrap_libsbml("document.createModel()")

        # find related objects
        # validate objects
        # group objects by model class
        grouped_objects = {}
        for obj in objects:
            obj_class = obj.__class__
            if obj_class not in grouped_objects:
                grouped_objects[obj_class] = []
            if obj not in grouped_objects[obj_class]:
                grouped_objects[obj_class].append(obj)

        # Dependencies among entities force an order on their creation.
        # E.g. compartments must be defined before the species that they contain
        model_order = [Submodel, Compartment, Species, Reaction]

        # add objects into SBMLDocument
        for model in model_order:
            if model in grouped_objects:
                for obj in grouped_objects[model]:
                    obj.add_to_sbml_doc(sbml_model)

        return document

    @staticmethod
    def read(document):
        """ Read a model in a libSBML `SBMLDocument` into a `wc_lang` model.

        Warning: `wc_lang` and SBML semantics are not equivalent. Thus, this
        SBML information will not be read:
            * xxx

        Args:
            document (:obj:`SBMLDocument`): representation of an SBML model

        Returns:
            :obj:`dict`: model objects grouped by `Model`

        Raises:
            :obj:`ValueError`: if ...
        """
        # TODO: write
        pass
