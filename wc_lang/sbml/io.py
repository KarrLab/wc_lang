""" Reading/writing a WC model to/from an SBML representation

Representations include
* Files
* Strings

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-09-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import libsbml
import os
import warnings
from os.path import split, join
from six import iteritems

from obj_model import Validator
from wc_lang.sbml.util import (init_sbml_model, create_sbml_doc_w_fbc)
import wc_lang

'''
wc_lang to SBML mapping to support FBA modeling
Individual wc_lang submodels that use dFBA are mapped to individual SBML documents and files.

WC                              SBML                                                                  Status
-----                           -----                                                                 ------
Model                           Model                                                                 Ignored
Taxon                           None, perhaps make SBML annotations                                   Ignored
Submodel                        Model                                                                 Implemented
DfbaObjective                   Objective                                                             Mostly Implemented
Compartment                     Compartment                                                           Implemented
SpeciesType                     SpeciesType aren't defined                                            NA
Species                         Species                                                               Implemented
DistributionInitConcentration   Distributions of initial concentrations are incorporated in Species   NA
Reaction                        Reaction, with FbcReactionPlugin for DFBA submodels                   Implemented
SpeciesCoefficient              SpeciesReference in a Reaction                                        Implemented
RateLaw                         KineticLaw                                                            Ignored
RateLawExpression
DfbaNetSpecies
DfbaNetReaction                 TBD
Parameter                       Parameter                                                             Implemented
Reference
DatabaseReference

wc_lang attribute to SBML mapping:

WC Model                SBML Model
--------                ----------
comments                notes
references              notes, as a Python dict
'''

#   class Reader(object):
#       """ Read model objects from an SBML representation """
#
#       @staticmethod
#       def run(self, objects, models, path=None, get_related=True,
#           title=None, description=None, keywords=None, version=None,
#           language=None, creator=None):
#           """ Write a list of model objects to a string or a file in an SBML format.
#
#           Warning: this is an unimplemented placeholder.
#
#           Args:
#               objects (:obj:`list`): list of objects
#               models (:obj:`list`): list of model, in the order that they should
#                   appear as worksheets; all models which are not in `models` will
#                   follow in alphabetical order
#               path (:obj:`str`, optional): path of SBML file to write
#               title (:obj:`str`, optional): title
#               description (:obj:`str`, optional): description
#               keywords (:obj:`str`, optional): keywords
#               version (:obj:`str`, optional): version
#               language (:obj:`str`, optional): language
#               creator (:obj:`str`, optional): creator
#           """
#           pass


class Writer(object):
    """ Write an SBML representation of a model  """

    @staticmethod
    def run(model, algorithms=None, path=None):
        """ Write the `model`'s submodels in SBML.

        Each `Submodel` in `Model` `model` whose algorithm is in `algorithms`
        is converted into a separate SBML document.
        If `path` is None, then the SBML is returned in string(s), otherwise it's written to file(s)
        which are named `path + submodel.id + suffix`.

        Args:
            model (:obj:`Model`): a `Model`
            algorithms (:obj:`list`, optional): list of `SubmodelAlgorithm` attributes, defaulting
                to `[SubmodelAlgorithm.dfba]`
            path (:obj:`str`, optional): prefix of path of SBML file(s) to write

        Returns:
            :obj:`dict` of `str`:
                if `path` is None, a dictionary SBML documents as strings, indexed by submodel ids,
                otherwise a list of SBML file(s) created
        """
        if algorithms is None:
            algorithms = [wc_lang.SubmodelAlgorithm.dfba]
        sbml_documents = {}
        for submodel in model.get_submodels():
            if submodel.algorithm in algorithms:
                objects = [submodel] + \
                    submodel.dfba_net_reactions + \
                    model.get_compartments() + \
                    submodel.get_species() + \
                    submodel.get_parameters() + \
                    submodel.reactions
                if submodel.dfba_obj:
                    objects.append(submodel.dfba_obj)
                sbml_documents[submodel.id] = SBMLExchange.write(objects)
        if not sbml_documents:
            raise ValueError("No submodel.algorithm in algorithms '{}'.".format(algorithms))
        if path is None:
            return sbml_documents
        else:
            ext = '.sbml'
            (dirname, basename) = split(path)
            files = []
            if not os.access(dirname, os.W_OK):
                raise ValueError("Writer.run() cannot write to directory '{}'.".format(dirname))
            for id, sbml_doc in iteritems(sbml_documents):
                dest = join(dirname, basename + '-' + id + ext)
                dest = str(dest)
                files.append(dest)
                if not libsbml.writeSBMLToFile(sbml_doc, dest):
                    raise ValueError("SBML document for submodel '{}' could not be written to '{}'.".format(
                        id, dest))
            return files


class SBMLExchange(object):
    """ Exchange `wc_lang` model to/from a libSBML SBML representation """

    @staticmethod
    def write(objects):
        """ Write the `wc_lang` model described by `objects` to a libSBML `SBMLDocument`.

        Warning: `wc_lang` and SBML semantics are not equivalent.

        Algorithm:
            * validate objects
            * (do not add related objects, as only certain model types can be written to SBML)
            * group objects by model class
            * add objects to the SBML document in dependent order

        Args:
            objects (:obj:`list`): list of objects

        Returns:
            :obj:`SBMLDocument`: an SBMLDocument containing `objects`

        Raises:
            :obj:`LibSBMLError`: if the SBMLDocument cannot be created
        """
        # Create an empty SBMLDocument object.
        sbml_document = create_sbml_doc_w_fbc()

        # Create the SBML Model object inside the SBMLDocument object.
        sbml_model = init_sbml_model(sbml_document)

        error = Validator().run(objects)
        if error:
            warnings.warn('Some data will not be written because objects are not valid:\n  {}'.format(
                str(error).replace('\n', '\n  ').rstrip()), UserWarning)

        grouped_objects = {}
        for obj in objects:
            obj_class = obj.__class__
            if obj_class not in grouped_objects:
                grouped_objects[obj_class] = []
            if obj not in grouped_objects[obj_class]:
                grouped_objects[obj_class].append(obj)

        # dependencies among libSBML model classes constrain the order
        # in which wc_lang classes must be written to SBML:
        #     Submodel depends on nothing
        #     Compartment depends on nothing
        #     Parameter depends on nothing
        #     Compartment must precede Species
        #     Compartment must precede Reaction
        #     Species must precede Reaction
        #     Species must precede DfbaNetReaction
        #     Reaction must precede DfbaObjective
        #     DfbaNetReaction must precede DfbaObjective
        # This partial order is satisfied by this sequence:
        model_order = [
            wc_lang.Submodel,
            wc_lang.Compartment,
            wc_lang.Parameter,
            wc_lang.Species,
            wc_lang.Reaction,
            wc_lang.DfbaNetReaction,
            wc_lang.DfbaObjective,
        ]

        # add objects into libsbml.SBMLDocument
        for model in model_order:
            if model in grouped_objects:
                for obj in grouped_objects[model]:
                    obj.add_to_sbml_doc(sbml_document)

        return sbml_document

    @staticmethod
    def write_submodel(submodel):
        """ Create a libSBML `SBMLDocument` containing `submodel`.

        To enable use of cobrapy to solve dFBA submodels, and avoid cumbersome SBML/libSBML
        submodels export one `wc_lang.Submodel` into one libSBML model.

        Args:
            submodel (:obj:`Submodel`): a submodel

        Returns:
            :obj:`libsbml.SBMLDocument`: an SBMLDocument containing `submodel` as a libSBML model

        Raises:
            :obj:`ValueError`: if the SBMLDocument cannot be created
        """
        objects = [submodel] + \
            submodel.model.get_compartments() + \
            submodel.get_species() + \
            submodel.reactions + \
            submodel.dfba_net_reactions + \
            submodel.model.get_parameters()
        if submodel.dfba_obj:
            objects.append(submodel.dfba_obj)

        return SBMLExchange.write(objects)

#       @staticmethod
#       def read(document):
#           """ Read a model in a libSBML `SBMLDocument` into a `wc_lang` model.
#
#           Warning: this is an unimplemented placeholder.
#           """
#           pass
