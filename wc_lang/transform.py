""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from operator import attrgetter
from wc_lang.core import Submodel
import itertools


class MergeAlgorithmicallyLikeSubmodels(object):
    """ Constructs models in which algorithmically-like submodels have been merged. """

    @staticmethod
    def transform(model):
        """ Construct a model in which algorithmically-like submodels have been merged.

        Args:
            model (:obj:`wc_lang.core.Model`): model definition

        Returns:
            :obj:`wc_lang.core.Model`: model with submodels of the same simulation algorithm merged
        """

        # copy model
        merged_model = model.copy()

        # group submodels by algorithms
        sorted_submodels = list(merged_model.submodels)
        sorted_submodels.sort(key=attrgetter('algorithm'))
        grouped_submodels = itertools.groupby(sorted_submodels, attrgetter('algorithm'))

        merged_submodels = set()
        for algorithm, group in grouped_submodels:
            submodels = tuple(group)

            # calculate id, name
            id = "_".join([submodel.id for submodel in submodels])
            name = "-".join([submodel.name for submodel in submodels])

            # instantiate merged submodel
            merged_submodel = Submodel(model=merged_model, id=id, name=name, algorithm=algorithm)

            # removed submodel from model; merge reactions, parameters, cross references, references
            for submodel in list(submodels):
                merged_model.submodels.remove(submodel)

                for rxn in list(submodel.reactions):
                    rxn.submodel = merged_submodel

                for param in list(submodel.parameters):
                    param.submodels.remove(submodel)
                    param.submodels.add(merged_submodel)

                for x_ref in list(submodel.cross_references):
                    x_ref.submodel = merged_submodel

                for ref in list(submodel.references):
                    ref.submodels.remove(submodel)
                    ref.submodels.add(merged_submodel)

        # return merged model
        return merged_model
