""" Encoding/decoding `wc_lang` models to/from SBML and
reading/writing SBML-encoded models to/from XML files.

Note, `wc_lang` and SBML do not have equivalent semantics.
Consequently, (a) `wc_lang` models must be exported as multiple
SBML models, one for each submodel and (b) SBML models which
utilize features, such as events, which are not supported by
`wc_lang` cannot be imported.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2019-03-21
:Copyright: 2017-2019, Karr Lab
:License: MIT
"""

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

from obj_model.utils import group_objects_by_model
from wc_lang.transform.prep_for_sbml import PrepForSbmlTransform
from wc_lang.sbml.util import LibSbmlInterface
from wc_utils.util.ontology import wcm_ontology, are_terms_equivalent
from wc_utils.util.units import unit_registry
import abc
import glob
import libsbml
import obj_model
import os
import warnings
import wc_lang.core


class SbmlWriter(object):
    """ Write `wc_lang` models to collections of SBML-encoded XML files, one for each submodel """

    def run(self, model, dirname):
        """ Write the submodels of a `wc_lang` model to separate SBML-encoded XML files.

        Args:
            model (:obj:`wc_lang.core.Model`): `wc_lang` model
            dirname (:obj:`str`): path to directory to save SBML-encoded XML files for each submodel
        """
        # validate model
        model = PrepForSbmlTransform().run(model.copy())
        error = wc_lang.core.Validator().run(model)
        if error:
            warnings.warn('Model is invalid: ' + str(error), wc_lang.core.WcLangWarning)

        # split submodels into separate models
        core, submodels = model.submodels.gen_models()
        all_models = [core] + submodels
        all_models_ids = ['core'] + [m.submodels[0].id for m in submodels]

        # create a directory to save SBML-encoded XML documents for model
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

        # encode models in SBML and save to XML file
        for model, model_id in zip(all_models, all_models_ids):
            # encode models in SBML
            sbml_doc = SbmlExporter.run(model)

            # save SBML-encoded model to XML file
            path = os.path.join(dirname, model_id + '.xml')
            if not LibSbmlInterface.call_libsbml(libsbml.writeSBMLToFile, sbml_doc, path, returns_int=True):
                raise ValueError("Submodel '{}' could not be written to SBML at '{}'.".format(
                    model_id, path))


class SbmlReader(object):
    """ Read `wc_lang` models from SBML-encoded XML files, one for each submodel """

    def run(self, dirname):
        """ Read `wc_lang` models from SBML-encoded XML files, one for each submodel

        Args:
            dirname (:obj:`str`): path to directory that contains SBML-encoded submodels of a model

        Returns:
            model (:obj:`wc_lang.core.Model`): `wc_lang` model
        """
        merged_model = None
        sbml_reader = LibSbmlInterface.call_libsbml(libsbml.SBMLReader)

        paths = glob.glob(os.path.join(dirname, '*.xml'))
        core_path = os.path.join(dirname, 'core.xml')
        if core_path in paths:
            paths.remove(core_path)
            paths.insert(0, core_path)

        for path in paths:
            # read model from XML file
            sbml_doc = LibSbmlInterface.call_libsbml(sbml_reader.readSBMLFromFile, path)
            LibSbmlInterface.raise_if_error(sbml_doc, 'Model could not be read from {}'.format(path))

            # convert SBML-encoded model to wc_lang
            model = SbmlImporter.run(sbml_doc)

            # merge models
            if merged_model is None:
                merged_model = model
            else:
                merged_model.merge(model)

        # return merged model
        return merged_model


class SbmlExporter(object):
    """ Encode a `wc_lang` model with at most 1 submodel into SBML """

    @classmethod
    def run(cls, model):
        """ Encode a `wc_lang` model with at most 1 submodel into SBML

        * Validate model
        * Create SBML document
        * Create SBML model
        * Encode model objects in SBML and add to SBML model in dependent order

        Args:
            model (:obj:`wc_lang.core.Model`): `wc_lang` model with at most 1 submodel

        Returns:
            :obj:`libsbml.Model`: SBML-encoded model
        """
        # verify model has at most 1 submodel
        if len(model.submodels) > 1:
            raise ValueError('Only 1 submodel can be encoded to SBML at a time')

        # validate model
        error = wc_lang.core.Validator().run(model)
        if error:
            warnings.warn('Model is invalid: ' + str(error), wc_lang.core.WcLangWarning)

        # determine SBML packages needed to export model
        packages = {}
        for submodel in model.submodels:
            if submodel.framework == wcm_ontology['WCM:dynamic_flux_balance_analysis']:
                packages['fbc'] = 2

        # Create an SBML document
        sbml_doc = LibSbmlInterface.create_doc(packages=packages)

        # Create a SBML model
        sbml_model = LibSbmlInterface.init_model(model, sbml_doc, packages=packages)

        # get model objects, grouped by type
        model_objs = group_objects_by_model(model.get_related())

        # add objects to SBML model
        # dependencies among libSBML model classes constrain the order
        model_order = []
        if model.submodels:
            model_order.append(wc_lang.core.Submodel)
        else:
            model_order.append(wc_lang.core.Model)
        model_order.extend([
            wc_lang.core.Compartment,
            wc_lang.core.Parameter,
            wc_lang.core.Species,
            wc_lang.core.Observable,
            wc_lang.core.Function,
            wc_lang.core.Reaction,
            wc_lang.core.DfbaObjReaction,
            wc_lang.core.DfbaObjective,
        ])
        for model in model_order:
            for obj in model_objs.get(model, []):
                obj.export_to_sbml(sbml_model)

        # verify document is compatible with SBML, valid SBML, and consistent
        LibSbmlInterface.verify_doc(sbml_doc)

        # return SBML document
        return sbml_doc


class SbmlImporter(object):
    """ Import a `wc_lang` model from an SBML-encoded model """

    @classmethod
    def run(cls, sbml_doc):
        """ Import a `wc_lang` model from an SBML-encoded model

        Args:
            sbml_doc (:obj:`libsbml.SBMLDocument`): SBML document with SBML-encoded model

        Returns:
            :obj:`wc_lang.core.Model`: `wc_lang` model
        """
        LibSbmlInterface.verify_doc(sbml_doc)
        sbml_model = LibSbmlInterface.call_libsbml(sbml_doc.getModel)

        # initialize model and submodel
        model = wc_lang.core.Model(created=None, updated=None)

        if LibSbmlInterface.call_libsbml(sbml_model.getNumReactions, returns_int=True):
            submodel = model.submodels.create()
            submodel.import_from_sbml(sbml_model)
        else:
            model.import_from_sbml(sbml_model)

        types = [
            ((wc_lang.core.Compartment, 'compartments'), (sbml_model, 'getNumCompartments', 'getCompartment')),
            ((wc_lang.core.Species, 'species'), (sbml_model, 'getNumSpecies', 'getSpecies')),
            ((wc_lang.core.Parameter, 'parameters'), (sbml_model, 'getNumParameters', 'getParameter')),
            ((wc_lang.core.Observable, 'observables'), (sbml_model, 'getNumRules', 'getRule')),
            ((wc_lang.core.Function, 'functions'), (sbml_model, 'getNumRules', 'getRule')),
            ((wc_lang.core.Reaction, 'reactions'), (sbml_model, 'getNumReactions', 'getReaction')),
        ]

        if LibSbmlInterface.call_libsbml(sbml_doc.isSetPackageRequired, 'fbc'):
            sbml_plugin = LibSbmlInterface.call_libsbml(sbml_model.getPlugin, 'fbc')
            types.append(((wc_lang.core.DfbaObjReaction, 'dfba_obj_reactions'), (sbml_model, 'getNumReactions', 'getReaction')))
            types.append(((wc_lang.core.DfbaObjective, 'dfba_objs'), (sbml_plugin, 'getNumObjectives', 'getObjective')))

        # create model objects
        idx_to_wc_lang_obj = {}
        id_to_wc_lang_obj = {}
        for (wc_lang_type, wc_lang_type_container), (sbml_root, sbml_type_counter, sbml_type_getter) in types:
            n_sbml_objs = LibSbmlInterface.call_libsbml(getattr(sbml_root, sbml_type_counter), returns_int=True)
            idx_to_wc_lang_obj[wc_lang_type] = {}
            id_to_wc_lang_obj[wc_lang_type] = {}
            for i_sbml_obj in range(n_sbml_objs):
                sbml_obj = LibSbmlInterface.call_libsbml(getattr(sbml_root, sbml_type_getter), i_sbml_obj)
                if not LibSbmlInterface.call_libsbml(sbml_obj.getIdAttribute).startswith(wc_lang_type.__name__ + '__'):
                    continue
                wc_lang_obj = getattr(model, wc_lang_type_container).create()
                wc_lang_obj.import_from_sbml(sbml_obj)
                idx_to_wc_lang_obj[wc_lang_type][i_sbml_obj] = wc_lang_obj
                id_to_wc_lang_obj[wc_lang_type][wc_lang_obj.id] = wc_lang_obj

        # link model objects
        if LibSbmlInterface.call_libsbml(sbml_model.getNumReactions, returns_int=True):
            submodel.import_relations_from_sbml(sbml_model, id_to_wc_lang_obj)
        else:
            model.import_relations_from_sbml(sbml_model, id_to_wc_lang_obj)
        for (wc_lang_type, wc_lang_type_container), (sbml_root, sbml_type_counter, sbml_type_getter) in types:
            for i_sbml_obj, wc_lang_obj in idx_to_wc_lang_obj[wc_lang_type].items():
                sbml_obj = getattr(sbml_root, sbml_type_getter)(i_sbml_obj)
                wc_lang_obj.import_relations_from_sbml(sbml_obj, id_to_wc_lang_obj)

        # return model
        return model
