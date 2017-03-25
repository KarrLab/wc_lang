""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from abc import ABCMeta, abstractmethod
from six import with_metaclass
from wc_lang.core import Model, Submodel, SubmodelAlgorithm, RateLawDirection
import itertools
import sys


def get_transforms():
    """ Get dictionary of available transform classes

    Returns:
        :obj:`dict` of `str`: `class`: dictionary of available transform classes
    """
    module = sys.modules[__name__]
    transforms = {}
    for attr_name in dir(module):
        attr = getattr(module, attr_name)
        if isinstance(attr, type) and issubclass(attr, Transform) and attr is not Transform:
            transforms[attr.Meta.id] = attr

    return transforms


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
        key_func = lambda submodel: submodel.algorithm.value
        sorted_submodels = sorted(model.submodels, key=key_func)
        grouped_submodels = itertools.groupby(sorted_submodels, key_func)

        for algorithm, group in grouped_submodels:
            submodels = tuple(group)

            # calculate id, name
            id = "_".join([submodel.id for submodel in submodels])
            name = "-".join([submodel.name for submodel in submodels])

            # instantiate merged submodel
            merged_submodel = Submodel(model=model, id=id, name=name, algorithm=SubmodelAlgorithm(algorithm))

            # removed submodel from model; merge reactions, parameters, cross references, references
            for submodel in submodels:
                model.submodels.remove(submodel)

                for rxn in submodel.reactions:
                    rxn.submodel = merged_submodel

                for param in reversed(submodel.parameters):
                    param.submodels.remove(submodel)
                    param.submodels.append(merged_submodel)

                for x_ref in reversed(submodel.cross_references):
                    x_ref.submodel = merged_submodel

                for ref in reversed(submodel.references):
                    ref.submodels.remove(submodel)
                    ref.submodels.append(merged_submodel)

        # return merged model
        return model


class SplitReversibleReactionsTransform(Transform):
    """ Split reversible reactions into separate forward and backward reactions """

    class Meta(object):
        id = 'SplitReversibleReactions'
        label = 'Split reversible reactions into separate forward and backward reactions'

    def run(self, model):
        """ Split reversible reactions into separate forward and backward reactions

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but with reversible reactions split into separate forward and backward reactions
        """

        for submodel in model.submodels:
            for rxn in submodel.reactions:
                if rxn.reversible:
                    # remove reversible reaction
                    submodel.reactions.remove(rxn)

                    # create separate forward and reverse reactions
                    rxn_for = submodel.reactions.create(
                        id='{}_forward'.format(rxn.id),
                        name='{} (forward)'.format(rxn.name),
                        reversible=False,
                        comments=rxn.comments,
                        references=rxn.references,
                    )
                    rxn_bck = submodel.reactions.create(
                        id='{}_backward'.format(rxn.id),
                        name='{} (backward)'.format(rxn.name),
                        reversible=False,
                        comments=rxn.comments,
                        references=rxn.references,
                    )

                    rxn.references = []

                    # copy participants and negate for backward reaction
                    for part in rxn.participants:
                        rxn_for.participants.append(part)
                        rxn_bck.participants.create(species=part.species, coefficient=-1 * part.coefficient)

                    rxn.participants = []

                    # copy rate laws
                    law_for = rxn.rate_laws.get(direction=RateLawDirection.forward)
                    law_bck = rxn.rate_laws.get(direction=RateLawDirection.backward)

                    if law_for:
                        law_for.reaction = rxn_for
                        law_for.direction = RateLawDirection.forward
                    if law_bck:
                        law_bck.reaction = rxn_bck
                        law_bck.direction = RateLawDirection.forward

                    # cross references
                    for x_ref in rxn.cross_references:
                        rxn_for.cross_references.create(
                            database=x_ref.database,
                            id=x_ref.id,
                            url=x_ref.url)

                        rxn_bck.cross_references.create(
                            database=x_ref.database,
                            id=x_ref.id,
                            url=x_ref.url)

                        rxn.cross_references.remove(x_ref)

        return model
