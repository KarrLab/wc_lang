""" Encoding/decoding `wc_lang` models to/from SBML and
reading/writing SBML-encoded models to/from XML files, one per submodel.

This enables WC-Lang-encoded models to be "round-tripped" with this sequence of operations:

1. Model is transformed to remove metadata that can't be exported to SBML
2. Model is split into multiple models, each with one submodel
3. Each model with one submodel is exported to SBML
4. Each model with one submodel is imported from SBML
5. Models with individual submodels are merged into a single model

WC-Lang models are exported to SBML with the following class mapping:

=============================  =====================
WC-Lang                        SBML
=============================  =====================
Model                          Model
Taxon                          Model.annotation
Environment                    Model.annotation
Submodel                       Model
Compartment                    Compartment
SpeciesType                    Species.annotation
Species                        Species
DistributionInitConcentration  Species.initialAmount
Observable                     AssignmentRule
ObservableExpression           AssignmentRule.math
Function                       AssignmentRule
FunctionExpression             AssignmentRule.math
Reaction                       Reaction
SpeciesCoefficient             SpeciesReference
RateLaw                        KineticLaw
RateLawExpression              KineticLaw.math
Parameter                      Parameter
DfbaObjective                  Objective
DfbaObjectiveExpression        FluxObjective
DfbaObjReaction                Reaction
DfbaObjSpecies                 SpeciesReference
StopCondition                  --
StopConditionExpression        --
Observation                       --
Conclusion                 --
Reference                      --
Author                         Model.annotation
Change                         --
Identifier              SBase.annotation
pint.Unit                      UnitDefinition
obj_model.Model.comments       SBase.notes
=============================  =====================

In addition, WC-Lang attributes which have no equivalent SBML attribute are mapped to
``SBase.annotation``.

Note, that because WC-Lang and SBML have different semantics, some of aspects of
WC-Lang cannot be mapped to SBML and vice-versa.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2019-03-21
:Copyright: 2017-2019, Karr Lab
:License: MIT
"""

from obj_model.utils import group_objects_by_model
from wc_lang.transform.prep_for_sbml import PrepForSbmlTransform
from wc_lang.sbml.util import LibSbmlInterface
from wc_onto import onto
from wc_utils.util.ontology import are_terms_equivalent
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

        Raises:
            :obj:`ValueError`: if the model could not be written to a SBML-encoded file
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

            import obj_model.core
            error = obj_model.core.Validator().run(model, get_related=True)
            assert error is None, str(error)


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
            :obj:`libsbml.SBMLDocument`: SBML document with SBML-encoded model

        Raises:
            :obj:`ValueError`: if the model cannot be exported to SBML because it contains multiple submodels
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
            if submodel.framework == onto['WC:dynamic_flux_balance_analysis']:
                packages['fbc'] = 2

        # create an SBML document
        sbml_doc = LibSbmlInterface.create_doc(packages=packages)

        # create a SBML model
        sbml_model = LibSbmlInterface.init_model(model, sbml_doc, packages=packages)

        # add objects to SBML model
        related_objs = model.get_related()
        sbml_objs = []
        for obj in related_objs:
            sbml_obj = obj.export_to_sbml(sbml_model)
            sbml_objs.append(sbml_obj)

        # add object relationships to SBML
        for obj, sbml_obj in zip(related_objs, sbml_objs):
            obj.export_relations_to_sbml(sbml_model, sbml_obj)

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
