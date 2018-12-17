""" Create implicit exchange reactions for dFBA submodels.

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang.core import SubmodelAlgorithm, ReactionFluxBoundUnit
import wc_lang.config.core


class CreateImplicitDfbaExchangeReactionsTransform(Transform):
    """ Create implicit exchange reactions for dFBA submodels.

    To enable FBA to represent a closed system, create implicit forward exchange reactions

    * Metabolites transported from the extracellular environment (generate a
      "-> species_type[e]" reaction for each extracellular species involved in each dFBA submodel)
    * TODO: Metabolites recycled from the byproducts of other submodels
    """

    class Meta(object):
        id = 'CreateImplicitDfbaExchangeReactions'
        label = 'Create implicit exchange reactions for dFBA submodels'

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        config = wc_lang.config.core.get_config()['wc_lang']
        ext_comp = model.compartments.get_one(id=config['EXTRACELLULAR_COMPARTMENT_ID'])
        rxn_id_template = config['dfba']['exchange_reaction_id_template']
        rxn_name_template = config['dfba']['exchange_reaction_name_template']
        ex_flux_bound_carbon = config['dfba']['ex_flux_bound_carbon']
        ex_flux_bound_no_carbon = config['dfba']['ex_flux_bound_no_carbon']

        for submodel in model.submodels:
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                for species in submodel.get_species():
                    if species.compartment == ext_comp:

                        rxn = submodel.reactions.create(
                            model=model,
                            id=rxn_id_template.format(submodel.id,
                                                      species.species_type.id,
                                                      species.compartment.id),
                            name=rxn_name_template.format(submodel.name,
                                                          species.species_type.name,
                                                          species.compartment.name),
                            reversible=True)

                        part = species.species_coefficients.get_one(coefficient=1.)
                        if part:
                            rxn.participants.append(part)
                        else:
                            rxn.participants.create(species=species, coefficient=1.)

                        if species.species_type.has_carbon():
                            rxn.flux_min = -ex_flux_bound_carbon
                            rxn.flux_max = ex_flux_bound_carbon
                        else:
                            rxn.flux_min = -ex_flux_bound_no_carbon
                            rxn.flux_max = ex_flux_bound_no_carbon
                        rxn.flux_bound_units = ReactionFluxBoundUnit['M s^-1']

        return model
