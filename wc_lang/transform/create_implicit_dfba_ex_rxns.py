""" Create implicit exchange reactions for dFBA submodels.

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_onto import onto
from wc_utils.util.ontology import are_terms_equivalent
from wc_utils.util.units import unit_registry
import wc_lang.config.core
import wc_lang.core


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
        ex_carbon = config['dfba']['flux_bounds']['ex_carbon']
        ex_no_carbon = config['dfba']['flux_bounds']['ex_no_carbon']

        carbon_flux_bounds = wc_lang.core.FluxBounds()
        carbon_flux_bounds.min = -ex_carbon
        carbon_flux_bounds.max = ex_carbon
        carbon_flux_bounds.units = unit_registry.parse_units('M s^-1')

        no_carbon_flux_bounds = wc_lang.core.FluxBounds()
        no_carbon_flux_bounds.min = -ex_no_carbon
        no_carbon_flux_bounds.max = ex_no_carbon
        no_carbon_flux_bounds.units = unit_registry.parse_units('M s^-1')

        for submodel in model.submodels:
            if are_terms_equivalent(submodel.framework, onto['WC:dynamic_flux_balance_analysis']):
                for species in submodel.get_children(kind='submodel', __type=wc_lang.core.Species):
                    if species.compartment == ext_comp:
                        id = rxn_id_template.format(submodel.id,
                                                    species.species_type.id,
                                                    species.compartment.id)
                        name = rxn_name_template.format(submodel.name or submodel.id,
                                                        species.species_type.name or species.species_type.id,
                                                        species.compartment.name or species.compartment.id)
                        participants = [species.species_coefficients.get_or_create(coefficient=1.)]
                        reversible = True
                        if species.species_type.structure and species.species_type.structure.has_carbon():
                            flux_bounds = carbon_flux_bounds
                        else:
                            flux_bounds = no_carbon_flux_bounds

                        rxn = model.reactions.get_one(id=id)
                        if rxn:
                            assert rxn.submodel == submodel
                            assert rxn.name == name
                            assert rxn.participants == participants
                            assert rxn.reversible == reversible
                            assert wc_lang.core.FluxBounds.Meta.attributes['min'].value_equal(
                                rxn.flux_bounds.min, flux_bounds.min)
                            assert wc_lang.core.FluxBounds.Meta.attributes['max'].value_equal(
                                rxn.flux_bounds.max, flux_bounds.max)
                            assert wc_lang.core.FluxBounds.Meta.attributes['units'].value_equal(
                                rxn.flux_bounds.units, flux_bounds.units)
                        else:
                            rxn = model.reactions.create(id=id)
                            rxn.submodel = submodel
                            rxn.name = name
                            rxn.participants = participants
                            rxn.reversible = reversible
                        rxn.flux_bounds = flux_bounds

        return model
