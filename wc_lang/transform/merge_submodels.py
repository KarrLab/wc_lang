""" Merge groups of algorithmically-like submodels into individual submodels.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang.core import (Model, Submodel, Reaction,
                          DfbaObjective, DfbaObjectiveExpression, DfbaObjReaction,
                          Observation, Evidence, Conclusion, Identifier, Reference,
                          Change)
from wc_onto import onto
import copy
import itertools


class MergeAlgorithmicallyLikeSubmodelsTransform(Transform):
    """ Merge groups of algorithmically-like submodels into individual submodels """

    class Meta(object):
        id = 'MergeAlgorithmicallyLikeSubmodels'
        label = 'Merge groups of algorithmically-like submodels into individual submodels'

    def run(self, model):
        """ Merge groups of algorithmically-like submodels into individual submodels

        * dFBA objectives are merged by summing

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but with submodels of the same simulation framework merged
        """

        # group submodels by frameworks
        def key_func(submodel):
            if submodel.dfba_obj:
                dfba_obj_units = submodel.dfba_obj.units
            else:
                dfba_obj_units = None
            return (submodel.framework.id, dfba_obj_units)
        sorted_submodels = sorted(model.submodels, key=key_func)
        grouped_submodels = itertools.groupby(sorted_submodels, key_func)

        for (framework, dfba_obj_units), group in grouped_submodels:
            submodels = tuple(group)

            # calculate id, name
            id = "_".join([submodel.id for submodel in submodels])
            name = "-".join([submodel.name for submodel in submodels])

            # instantiate merged submodel
            merged_submodel = Submodel(model=model, id=id, name=name, framework=onto[framework])

            if  framework == 'WC:dynamic_flux_balance_analysis':
                merged_dfba_obj = merged_submodel.dfba_obj = model.dfba_objs.create(
                    name='dFBA objective ({})'.format(', '.join(submodel.name for submodel in submodels)),
                    units=dfba_obj_units)
                merged_dfba_obj.id = merged_dfba_obj.gen_id()

                merged_dfba_expression = []
                objs_for_merged_dfba_expression = {
                    Reaction: {},
                    DfbaObjReaction: {},
                }

            # removed submodel from model
            # merge submodels
            # - model
            # - identifiers
            # - evidence
            # - conclusions
            # - references
            # - reactions
            # - dfba_obj
            # - dfba_obj_reactions
            for submodel in submodels:
                # assert that all types of related objects will be merged
                assert set(attr.related_class for attr in Submodel.Meta.local_attributes.values() if attr.related_class) == set(
                    [Model, Evidence, Conclusion, Identifier, Reference, Reaction,
                     DfbaObjective, DfbaObjReaction, Change])

                model.submodels.remove(submodel)

                for ev in list(submodel.evidence):
                    submodel.evidence.remove(ev)
                    merged_submodel.evidence.append(ev)

                for conclusion in list(submodel.conclusions):
                    submodel.conclusions.remove(conclusion)
                    merged_submodel.conclusions.append(conclusion)

                for identifier in list(submodel.identifiers):
                    submodel.identifiers.remove(identifier)
                    merged_submodel.identifiers.append(identifier)

                for ref in list(submodel.references):
                    submodel.references.remove(ref)
                    merged_submodel.references.append(ref)

                for rxn in list(submodel.reactions):
                    submodel.reactions.remove(rxn)
                    merged_submodel.reactions.append(rxn)

                if submodel.dfba_obj:
                    # assert that all types of related objects will be merged
                    assert set(attr.related_class for attr in DfbaObjective.Meta.local_attributes.values() if attr.related_class) == set(
                        [Model, Submodel, Evidence, Conclusion, Identifier, Reference, DfbaObjectiveExpression])

                    model.dfba_objs.remove(submodel.dfba_obj)

                    for ev in list(submodel.dfba_obj.evidence):
                        submodel.dfba_obj.evidence.remove(ev)
                        merged_submodel.dfba_obj.evidence.append(ev)

                    for conclusion in list(submodel.dfba_obj.conclusions):
                        submodel.dfba_obj.conclusions.remove(conclusion)
                        merged_submodel.dfba_obj.conclusions.append(conclusion)

                    for identifier in list(submodel.dfba_obj.identifiers):
                        submodel.dfba_obj.identifiers.remove(identifier)
                        merged_submodel.dfba_obj.identifiers.append(identifier)

                    for ref in list(submodel.dfba_obj.references):
                        submodel.dfba_obj.references.remove(ref)
                        merged_submodel.dfba_obj.references.append(ref)

                    if submodel.dfba_obj.expression:
                        assert set(attr.related_class for attr in DfbaObjectiveExpression.Meta.local_attributes.values()
                                   if attr.related_class) == set([DfbaObjective, Reaction, DfbaObjReaction])

                        if submodel.dfba_obj.expression.expression:
                            merged_dfba_expression.append(submodel.dfba_obj.expression.expression)
                        for rxn in list(submodel.dfba_obj.expression.reactions):
                            submodel.dfba_obj.expression.reactions.remove(rxn)
                            objs_for_merged_dfba_expression[Reaction][rxn.id] = rxn
                        for dfba_obj_rxn in list(submodel.dfba_obj.expression.dfba_obj_reactions):
                            submodel.dfba_obj.expression.dfba_obj_reactions.remove(dfba_obj_rxn)
                            objs_for_merged_dfba_expression[DfbaObjReaction][dfba_obj_rxn.id] = dfba_obj_rxn

                for dfba_obj_rxn in list(submodel.dfba_obj_reactions):
                    submodel.dfba_obj_reactions.remove(dfba_obj_rxn)
                    merged_submodel.dfba_obj_reactions.append(dfba_obj_rxn)

            if framework == 'WC:dynamic_flux_balance_analysis':
                merged_dfba_obj.expression, error = DfbaObjectiveExpression.deserialize(
                    ' + '.join(merged_dfba_expression),
                    objs_for_merged_dfba_expression)
                assert error is None

        # return merged model
        return model
