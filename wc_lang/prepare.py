""" Prepare a model for simulation.

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from math import isnan
from obj_model.utils import get_component_by_id
from wc_lang.core import Submodel, SubmodelAlgorithm, Concentration, ConcentrationUnit
from wc_lang.transform.core import Transform
import wc_lang.config.core

# configuration
config = wc_lang.config.core.get_config()['wc_lang']


class PrepareModelTransform(Transform):
    """ Prepare a model for simulation by making implicit information in the model
    explicit.

    * Create explicit zero mean concentrations for species that don't have explicit
      mean concentrations (i.e. implicit zero concentrations)
    * Create implicit exchange reactions for dFBA submodels
    * Clip the flux bounds for the reactions in dFBA submodels to the default
      flux range
    """

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        CreateImplicitZeroConcentrationsTransform().run(model)
        CreateImplicitDfbaExchangeReactionsTransform().run(model)
        SetFiniteDfbaFluxBoundsTransform().run(model)


class CreateImplicitZeroConcentrationsTransform(Transform):
    """ Create the implicit zero concentrations in a model
    """

    class Meta(object):
        id = 'CreateImplicitZeroConcentrations'
        label = 'Create the implicit zero concentrations in a model'

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        for species in model.get_species():
            if species.concentration is None:
                species.concentration = Concentration(
                    id=Concentration.gen_id(species.id),
                    model=model,
                    species=species,
                    mean=0.0, std=0.0, units=ConcentrationUnit.molecules)

        return model


class CreateImplicitDfbaExchangeReactionsTransform(Transform):
    """ Create implicit exchange reactions for dFBA submodels.

        To enable FBA to represent a closed system, create implicit forward exchange reactions

        * Metabolites transported from the extracellular environment (generate a 
          "-> species_type[e]" reaction for each extracellular species involved in each dFBA submodel)
        * TODO: Metabolites recycled from the byproducts of other submodels
    """

    class Meta(object):
        id = 'CreateImplicitZeroConcentrations'
        label = 'Create implicit exchange reactions for dFBA submodels'

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        ext_comp = model.compartments.get_one(id=config['EXTRACELLULAR_COMPARTMENT_ID'])
        rxn_id_template = config['dfba']['exchange_reaction_id_template']
        rxn_name_template = config['dfba']['exchange_reaction_name_template']

        for submodel in model.submodels:
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                for species in submodel.get_species():
                    if species.compartment == ext_comp:
                        rxn = submodel.reactions.create(
                            id=rxn_id_template.format(submodel.id,
                                                      species.species_type.id,
                                                      species.compartment.id),
                            name=rxn_name_template.format(submodel.name,
                                                          species.species_type.name,
                                                          species.compartment.name),
                            model=model,
                            reversible=True)
                        rxn.participants.create(species=species, coefficient=1.)
        return model


class SetFiniteDfbaFluxBoundsTransform(Transform):
    """ Clip the flux bounds for the reactions in dFBA submodels to the 
    default flux range because some linear programming solvers require
    finite minimum and maximum flux bounds.

    The default bounds are provided in the wc_lang configuration.

    Specifically, the default minimum and maximum flux bounds are applied as follows:

        * reversible reactions:

          * min_flux = max(min_flux, min_reversible_flux_bound)
          * max_flux = min(max_flux, max_flux_bound)

        * irreversible reactions:

          * min_flux = max(min_flux, min_irreversible_flux_bound)
          * max_flux = min(max_flux, max_flux_bound)
    """

    class Meta(object):
        id = 'SetFiniteDfbaFluxBoundsTransform'
        label = ('Clip the flux bounds for the reactions in dFBA submodels'
                 ' to the default flux range')

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        min_reversible_flux_bound = config['dfba']['min_reversible_flux_bound']
        min_irreversible_flux_bound = config['dfba']['min_irreversible_flux_bound']
        max_flux_bound = config['dfba']['max_flux_bound']

        for submodel in model.submodels:
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                for rxn in submodel.reactions:
                    if rxn.reversible:
                        min_flux = min_reversible_flux_bound
                    else:
                        min_flux = min_irreversible_flux_bound
                    max_flux = max_flux_bound

                    if rxn.min_flux is None or isnan(rxn.min_flux):
                        rxn.min_flux = min_flux
                    else:
                        rxn.min_flux = max(rxn.min_flux, min_flux)

                    if rxn.max_flux is None or isnan(rxn.max_flux):
                        rxn.max_flux = max_flux
                    else:
                        rxn.max_flux = min(rxn.max_flux, max_flux)
        return model
