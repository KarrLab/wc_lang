""" Utilities

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from obj_model.core import get_models as base_get_models
import wc_lang.core


def get_model_size(model):
    """ Get numbers of model components

    Args:
        model (:obj:`wc_lang.core.Model`): model

    Returns:
        :obj:`dict`: dictionary with numbers of each type of model component
    """
    return {
        "submodels": len(model.get_submodels()),
        "compartments": len(model.get_compartments()),
        "species_types": len(model.get_species_types()),
        "species": len(model.get_species()),
        "parameters": len(model.get_parameters()),
        "references": len(model.get_references()),
        "reactions": len(model.get_reactions()),
    }


def get_model_summary(model):
    """ Get textual summary of a model

    Args:
        model (:obj:`wc_lang.core.Model`): model

    Returns:
        :obj:`str`: textual summary of the model
    """
    return "Model with:" \
        + "\n{:d} submodels".format(len(model.get_submodels())) \
        + "\n{:d} compartments".format(len(model.get_compartments())) \
        + "\n{:d} species types".format(len(model.get_species_types())) \
        + "\n{:d} species".format(len(model.get_species())) \
        + "\n{:d} parameters".format(len(model.get_parameters())) \
        + "\n{:d} references".format(len(model.get_references())) \
        + "\n{:d} biomass reactions".format(len(model.get_biomass_reactions())) \
        + "\n{:d} reactions".format(len(model.get_reactions())) \
        + "\n{:d} rate laws".format(len(model.get_rate_laws()))


def get_reaction_string(reaction):
    """ Generate string representation of reaction stoichometry.

    Returns:
        :obj:`str`: string representation of reaction stoichometry
    """
    attr = reaction.__class__.Meta.attributes['participants']
    return attr.serialize(reaction.participants)


def get_models(inline=True):
    """ Get list of models
    Args:
        inline (:obj:`bool`, optional): if true, return inline models

    Returns:
        :obj:`list` of `class`: list of models
    """

    return base_get_models(module=wc_lang.core, inline=inline)
