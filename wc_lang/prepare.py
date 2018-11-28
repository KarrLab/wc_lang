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
import wc_lang.config.core

# configuration
config = wc_lang.config.core.get_config()['wc_lang']


class PrepareModel(object):
    """ Prepare a model for simulation by adding implicit information that is not explicitly
    specified in model definitions.

    * Insert 0 mean concentrations for species that don't have mean concentrations
    * Create implicit exchange reactions for dFBA submodels
    * Apply default flux bounds to the reactions in dFBA submodels
    """

    def __init__(self, model):
        self.model = model

    def run(self):
        """ Statically prepare a model for simulation by adding implicit information that is not explicitly
        specified in model definitions.
        """
        self.create_implicit_zero_concentrations()
        self.create_implicit_dfba_exchange_reactions()
        self.set_finite_dfba_flux_bounds()

    def create_implicit_zero_concentrations(self):
        """ Define implicit zero mean concentrations """
        model = self.model
        for species in model.get_species():
            if species.concentration is None:
                species.concentration = Concentration(
                    id=Concentration.gen_id(species.id),
                    model=model,
                    species=species,
                    mean=0.0, std=0.0, units=ConcentrationUnit.molecules)

    def create_implicit_dfba_exchange_reactions(self):
        """ Create implicit exchange reactions for dFBA submodels.

        To enable FBA to represent a closed system, create implicit forward exchange reactions

        * Metabolites transported from the extracellular environment (generate a 
          "-> species_type[e]" reaction for each extracellular species involved in each dFBA submodel)
        * TODO: Metabolites recycled from the byproducts of other submodels
        """
        model = self.model
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

    def set_finite_dfba_flux_bounds(self):
        """ Apply default flux bounds to the reactions in dFBA submodels because 
        some linear programming solvers require finite minimum and maximum flux bounds.

        The default default bounds are provided in the wc_lang configuration.

        Specifically, the default minimum and maximum flux bounds are applied as follows:

            * reversible reactions:

              * min_flux = max(min_flux, min_reversible_flux_bound)
              * max_flux = min(max_flux, max_flux_bound)

            * irreversible reactions:

              * min_flux = max(min_flux, min_irreversible_flux_bound)
              * max_flux = min(max_flux, max_flux_bound)
        """
        min_reversible_flux_bound = config['dfba']['min_reversible_flux_bound']
        min_irreversible_flux_bound = config['dfba']['min_irreversible_flux_bound']
        max_flux_bound = config['dfba']['max_flux_bound']

        model = self.model
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
