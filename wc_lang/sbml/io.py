""" Reading/writing a WC model to/from an SBML representation

Representations include
* Files
* Strings

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-24
:Copyright: 2017-2019, Karr Lab
:License: MIT
"""

from obj_model import Validator
from six import iteritems
from wc_lang.sbml.util import init_sbml_model, create_sbml_doc_w_fbc
from wc_utils.util.ontology import wcm_ontology
import libsbml
import os
import warnings
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
DfbaObjSpecies
DfbaObjReaction                 TBD
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


class SbmlWriter(object):
    """ Write an SBML representation of a model  """

    @classmethod
    def run(cls, model, out_dir):
        """ Write the submodels of a model to separate SBML-encoded XML files.

        Each :obj:`Submodel` in :obj:`Model` :obj:`model` is converted to a separate SBML document
        and saved to a separate XML file.

        Args:
            model (:obj:`Model`): a `Model`
            out_dir (:obj:`str`): path to directory to save SBML file(s)

        Returns:
            :obj:`dict`: dictionary that maps submodels to the paths where they are exported
        """
        # create a directory to save SBML-encoded XML documents for each submodel
        if not os.path.isdir(out_dir):
            os.mkdirs(out_dir)

        # save submodels to SBML-encoded XML files
        sbml_paths = {}
        for submodel in model.submodels:
            if submodel.framework != wcm_ontology['WCM:dynamic_flux_balance_analysis']:  # todo: remove
                continue
            sbml_paths[submodel] = cls.run_submodel(submodel, out_dir)
        return sbml_paths

    @classmethod
    def run_submodel(cls, submodel, out_dir):
        out_file = os.path.join(out_dir, submodel.id + '.xml')
        sbml_doc = SbmlConverter.run_submodel(submodel)
        if not libsbml.writeSBMLToFile(sbml_doc, out_file):
            raise ValueError("SBML document for submodel '{}' could not be written to '{}'.".format(
                submodel.id, out_file))
        return out_file


class SbmlConverter(object):
    """ Convert `wc_lang` model to/from a libSBML SBML representation """

    @classmethod
    def run(cls, model):
        sbml_docs = {}
        for submodel in model.submodels:
            if submodel.framework != wcm_ontology['WCM:dynamic_flux_balance_analysis']:  # todo: remove
                continue
            sbml_doc = cls.run_submodel(submodel)
            sbml_docs[submodel] = sbml_doc
        return sbml_docs

    @classmethod
    def run_submodel(cls, submodel):
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
        objects = [submodel] \
            + submodel.get_compartments() \
            + submodel.get_species() \
            + submodel.get_parameters() \
            + submodel.reactions
        if submodel.framework == wcm_ontology['WCM:dynamic_flux_balance_analysis']:
            objects.append(submodel.dfba_obj)
            objects.extend(submodel.dfba_obj_reactions)

        # DistributionInitConcentration
        # observables
        # functions

        #core.Taxon, core.Environment,
        #core.Evidence, core.Interpretation, core.Reference, core.Author, core.Change,

        return cls.write(objects)

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
            :obj:`LibSbmlError`: if the SBMLDocument cannot be created
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
        #     Species must precede DfbaObjReaction
        #     Reaction must precede DfbaObjective
        #     DfbaObjReaction must precede DfbaObjective
        # This partial order is satisfied by this sequence:
        model_order = [
            wc_lang.Submodel,
            wc_lang.Compartment,
            wc_lang.Parameter,
            wc_lang.Species,
            wc_lang.Reaction,
            wc_lang.DfbaObjReaction,
            wc_lang.DfbaObjective,
        ]

        # add objects into libsbml.SBMLDocument
        for model in model_order:
            if model in grouped_objects:
                for obj in grouped_objects[model]:
                    obj.add_to_sbml_doc(sbml_document)

        return sbml_document

#       @staticmethod
#       def read(document):
#           """ Read a model in a libSBML `SBMLDocument` into a `wc_lang` model.
#
#           Warning: this is an unimplemented placeholder.
#           """
#           pass
