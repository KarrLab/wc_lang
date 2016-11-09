""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2017-07-25
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang.core import Submodel
import copy
import itertools
import operator


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
        merged_model = copy.deepcopy(model)

        # group submodels by algorithms
        merged_submodels = []
        for algorithm, group in itertools.groupby(merged_model.submodels, operator.attrgetter('algorithm')):
            submodels = tuple(group)

            # calculate id, name
            id = "-".join([submodel.id for submodel in submodels])
            name = "-".join([submodel.name for submodel in submodels])

            # instantiate merged submodel
            merged_submodel = Submodel(id, name, algorithm, species=[], reactions=[], parameters=[])

            # merge species, reactions, parameters
            for submodel in submodels:
                for species in submodel.species:
                    merged_submodel.species.append(species)

                for reaction in submodel.reactions:
                    reaction.submodel = merged_model
                    merged_submodel.reactions.append(reaction)

                for parameter in submodel.parameters:
                    parameter.submodel = merged_submodel
                    merged_submodel.parameters.append(parameter)

            # get unique set of species
            unique_species = {}
            for species in merged_submodel.species:
                unique_species[species.id] = species
            merged_submodel.species = unique_species.values()

            # append to list of merged submodels
            merged_submodels.append(merged_submodel)

        # replace submodels with merged versions
        merged_model.submodels = merged_submodels

        # return merged model
        return merged_model
