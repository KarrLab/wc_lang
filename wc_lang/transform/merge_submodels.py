""" Merge groups of algorithmically-like submodels into individual submodels.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang import Model, Submodel, SubmodelAlgorithm
import copy
import itertools


class MergeAlgorithmicallyLikeSubmodelsTransform(Transform):
    """ Merge groups of algorithmically-like submodels into individual submodels """

    class Meta(object):
        id = 'MergeAlgorithmicallyLikeSubmodels'
        label = 'Merge groups of algorithmically-like submodels into individual submodels'

    def run(self, model):
        """ Merge groups of algorithmically-like submodels into individual submodels

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but with submodels of the same simulation algorithm merged
        """

        # group submodels by algorithms
        def key_func(submodel): return submodel.algorithm.value
        sorted_submodels = sorted(model.submodels, key=key_func)
        grouped_submodels = itertools.groupby(sorted_submodels, key_func)

        for algorithm, group in grouped_submodels:
            submodels = tuple(group)

            # calculate id, name
            id = "_".join([submodel.id for submodel in submodels])
            name = "-".join([submodel.name for submodel in submodels])

            # instantiate merged submodel
            merged_submodel = Submodel(model=model, id=id, name=name, algorithm=SubmodelAlgorithm(algorithm))

            # removed submodel from model; merge reactions, parameters, database references, references
            for submodel in submodels:
                model.submodels.remove(submodel)

                if submodel.dfba_obj:
                    dfba_obj = copy.copy(submodel.dfba_obj)
                    dfba_obj.submodel = merged_submodel

                for dfba_net_rxn in copy.copy(submodel.dfba_net_reactions):
                    dfba_net_rxn.submodel = merged_submodel

                for rxn in copy.copy(submodel.reactions):
                    rxn.submodel = merged_submodel

                for db_ref in copy.copy(submodel.db_refs):
                    db_ref.submodels.remove(submodel)
                    db_ref.submodels.append(merged_submodel)

                for ref in copy.copy(submodel.references):
                    ref.submodels.remove(submodel)
                    ref.submodels.append(merged_submodel)

        # return merged model
        return model
