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

from wc_lang.sbml.util import LibSbmlInterface
from wc_utils.util.ontology import wcm_ontology
import abc
import libsbml
import obj_model
import os
import warnings
import wc_lang.core

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


class SbmlWriter(object):
    """ Write `wc_lang` models to a a collection of SBML-encoded XML files, one for each submodel """

    @classmethod
    def run(cls, model, out_dir):
        """ Write the submodels of a `wc_lang` model to separate SBML-encoded XML files.

        Args:
            model (:obj:`wc_lang.core.Model`): `wc_lang` model
            out_dir (:obj:`str`): path to directory to save SBML-encoded XML files for each submodel

        Returns:
            :obj:`dict`: dictionary that maps submodels to the paths where they were exported
        """
        # encode submodels in SBML
        sbml_submodels = SbmlExporter.run(model)

        # create a directory to save SBML-encoded XML documents for each submodel
        if not os.path.isdir(out_dir):
            os.mkdirs(out_dir)

        # save submodels to SBML-encoded XML files
        sbml_submodel_paths = {}
        for submodel, sbml_submodel in sbml_submodels.items():
            sbml_submodel_paths[submodel] = cls.run_submodel(submodel, sbml_submodel, out_dir)
        return sbml_submodel_paths

    @classmethod
    def run_submodel(cls, submodel, sbml_submodel, out_dir):
        """ Write a `wc_lang` submodel to an SBML-encoded XML file

        Args:
            submodel (:obj:`wc_lang.core.Submodel`): `wc_lang` submodel
            sbml_submodel (:obj:`libsbml.SBMLDocument): SBML-encoded submodel
            out_dir (:obj:`str`): path to directory to save SBML-encoded XML file

        Returns:
            :obj:`str`: path to SBML-encoded XML file
        """
        out_file = os.path.join(out_dir, submodel.id + '.xml')
        if not libsbml.writeSBMLToFile(sbml_submodel, out_file):
            raise ValueError("SBML document for submodel '{}' could not be written to '{}'.".format(
                submodel.id, out_file))
        return out_file


class SbmlExporter(object):
    """ Export each submodel of a `wc_lang` model to an SBML-encoded model """

    @classmethod
    def run(cls, model):
        """ Export each submodel of a `wc_lang` model to an SBML-encoded model

        Args:
            model (:obj:`wc_lang.core.Model`): `wc_lang` model

        Returns:
            :obj:`dict`: dictionary that maps submodels to SBML-encoded copies of each submodel
        """

        # validate model
        error = wc_lang.core.Validator().run(model)
        if error:
            warnings.warn('Model is invalid: ' + str(error), UserWarning)

        # convert each submodel to an SBML model
        sbml_submodels = {}
        for submodel in model.submodels:
            if submodel.framework != wcm_ontology['WCM:dynamic_flux_balance_analysis']:
                continue
            sbml_submodels[submodel] = SubmodelSbmlExporter.run(submodel)
        return sbml_submodels


class SubmodelSbmlExporter(object):
    """ Export a `wc_lang` submodel to an SBML-encoded model """
    @classmethod
    def run(cls, submodel):
        """ Export a `wc_lang` submodel to a :obj:`libsbml.SBMLDocument`.        

        * Validate objects
        * (do not add related objects, as only certain model types can be written to SBML)
        * Group objects by model class
        * Add objects to the SBML document in dependent order

        Warning: `wc_lang` and SBML semantics are not equivalent.

        Args:
            submodel (obj:`wc_lang.core.Submodel`): `wc_lang` submodel

        Returns:
            :obj:`libsbml.SBMLDocument`: SBML-encoded submodel
        """
        packages = {}
        if submodel.framework == wcm_ontology['WCM:dynamic_flux_balance_analysis']:
            packages['fbc'] = 2

        # Create an empty libsbml.SBMLDocument object.
        sbml_doc = LibSbmlInterface.create_doc(packages=packages)

        # Create the SBML Model object inside the libsbml.SBMLDocument object.
        sbml_model = LibSbmlInterface.init_model(None, sbml_doc, packages=packages)

        objects = cls.get_submodel_objects(submodel)

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
            wc_lang.core.Submodel,
            wc_lang.core.Compartment,
            wc_lang.core.Parameter,
            wc_lang.core.Species,
            wc_lang.core.Reaction,
            wc_lang.core.DfbaObjReaction,
            wc_lang.core.DfbaObjective,
        ]

        # add objects into libsbml.SBMLDocument
        for model in model_order:
            if model in grouped_objects:
                for obj in grouped_objects[model]:
                    obj.add_to_sbml_model(sbml_model)

        return sbml_doc

    @classmethod
    def get_submodel_objects(cls, submodel):
        """ Get the children that must be exported with a submodel and validate that
        they are are a valid model

        Args:
            submodel (:obj:`wc_lang.core.Submodel`): `wc_lang` submodel

        Returns:
            :obj:`list` of :obj:`obj_model.Model`: `wc_lang` submodel and its `wc_lang` children
        """
        objects = [submodel] \
            + submodel.get_children(kind='submodel', __type=wc_lang.core.Compartment) \
            + submodel.get_children(kind='submodel', __type=wc_lang.core.Species) \
            + submodel.get_children(kind='submodel', __type=wc_lang.core.Parameter) \
            + submodel.reactions
        if submodel.framework == wcm_ontology['WCM:dynamic_flux_balance_analysis']:
            if submodel.dfba_obj:
                objects.append(submodel.dfba_obj)
            objects.extend(submodel.dfba_obj_reactions)

        # DistributionInitConcentration
        # observables
        # functions

        #core.Taxon, core.Environment,
        #core.Evidence, core.Interpretation, core.Reference, core.Author, core.Change,

        # validate objects
        error = obj_model.Validator().run(objects)
        if error:
            warnings.warn('Some data will not be written because objects are not valid:\n  {}'.format(
                str(error).replace('\n', '\n  ').rstrip()), UserWarning)

        # return objects
        return objects


class SbmlReader(object):
    """ Read `wc_lang` models from SBML-encoded XML files """

    @classmethod
    def run(self, out_dir,
            title=None, description=None, keywords=None, version=None,
            language=None, creator=None):
        """ Read `wc_lang` models from SBML-encoded XML files

        Args:
            out_dir (:obj:`str`): path to directory that contains SBML-encoded submodels of a model
            title (:obj:`str`, optional): title
            description (:obj:`str`, optional): description
            keywords (:obj:`str`, optional): keywords
            version (:obj:`str`, optional): version
            language (:obj:`str`, optional): language
            creator (:obj:`str`, optional): creator
        """
        pass


class SbmlImporter(object):
    """ Import a wc-lang-encoded model from an SBML-encoded model """

    @classmethod
    def run(cls, sbml_submodels):
        """ Import a `wc_lang` model from SBML-encoded submodels

        Args:
            sbml_submodels (:obj:`list` of :obj:`libsbml.SBMLDocument`): SBML-encoded submodels

        Returns:
            :obj:`wc_lang.core.Model`: `wc_lang` model
        """
        pass
