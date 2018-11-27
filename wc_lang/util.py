""" Utilities

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from obj_model import get_models as base_get_models
from wc_lang import core
from wc_utils.util import git


def get_model_size(model):
    """ Get numbers of model components

    Args:
        model (:obj:`core.Model`): model

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
        model (:obj:`core.Model`): model

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
        + "\n{:d} dFBA net reactions".format(len(model.get_dfba_net_reactions())) \
        + "\n{:d} reactions".format(len(model.get_reactions())) \
        + "\n{:d} rate laws".format(len(model.get_rate_laws()))


def get_models(inline=True):
    """ Get list of models
    Args:
        inline (:obj:`bool`, optional): if true, return inline models

    Returns:
        :obj:`list` of :obj:`class`: list of models
    """

    return base_get_models(module=core, inline=inline)


def set_git_repo_metadata_from_path(model, path='.'):
    """ Use Git to set the Git repository URL, branch, and revision metadata for a model

    Args:
        model (:obj:`core.Model`): model
        path (:obj:`str`, optional): path to the Git repository for the model
    """
    md = git.get_repo_metadata(dirname=path)
    model.url = md.url
    model.branch = md.branch
    model.revision = md.revision
