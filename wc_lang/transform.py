""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from abc import ABCMeta, abstractmethod
from operator import attrgetter
from six import with_metaclass
from wc_lang.core import Model, Submodel
import itertools


class Transform(with_metaclass(ABCMeta, object)):

    @abstractmethod
    def run(self, model):
        """ Transform a model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: transformed model
        """
        pass


class MergeAlgorithmicallyLikeSubmodelsTransform(Transform):
    """ Merge groups of algorithmically-like submodels into individual submodels """

    def run(self, model):
        """ Merge groups of algorithmically-like submodels into individual submodels

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but with submodels of the same simulation algorithm merged
        """

        # group submodels by algorithms
        sorted_submodels = list(model.submodels)
        sorted_submodels.sort(key=attrgetter('algorithm'))
        grouped_submodels = itertools.groupby(sorted_submodels, attrgetter('algorithm'))

        merged_submodels = set()
        for algorithm, group in grouped_submodels:
            submodels = tuple(group)

            # calculate id, name
            id = "_".join([submodel.id for submodel in submodels])
            name = "-".join([submodel.name for submodel in submodels])

            # instantiate merged submodel
            merged_submodel = Submodel(model=model, id=id, name=name, algorithm=algorithm)

            # removed submodel from model; merge reactions, parameters, cross references, references
            for submodel in list(submodels):
                model.submodels.remove(submodel)

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
        return model


class SplitReversibleReactionsTransform(Transform):
    """ Split reversible reactions into separate forward and backward reactions """

    def run(self, model):
        """ Split reversible reactions into separate forward and backward reactions

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but with reversible reactions split into separate forward and backward reactions
        """

        return model
