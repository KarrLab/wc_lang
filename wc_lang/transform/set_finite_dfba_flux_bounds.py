""" Clip the flux bounds for the reactions in dFBA submodels

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from .core import Transform
from math import isnan
from wc_onto import onto
from wc_utils.util.ontology import are_terms_equivalent
from wc_utils.util.units import unit_registry
import wc_lang.config.core
import wc_lang.core


class SetFiniteDfbaFluxBoundsTransform(Transform):
    """ Clip the flux bounds for the reactions in dFBA submodels to the
    default flux range because some linear programming solvers require
    finite minimum and maximum flux bounds.

    The default bounds are provided in the wc_lang configuration.

    Specifically, the default minimum and maximum flux bounds are applied as follows:

        * reversible reactions:

          * flux_min = max(flux_min, flux_min_bound_reversible)
          * flux_max = min(flux_max, flux_max_bound)

        * irreversible reactions:

          * flux_min = max(flux_min, flux_min_bound_irreversible)
          * flux_max = min(flux_max, flux_max_bound)
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
        config = wc_lang.config.core.get_config()['wc_lang']
        flux_min_bound_reversible = config['dfba']['flux_bounds']['min_reversible']
        flux_min_bound_irreversible = config['dfba']['flux_bounds']['min_irreversible']
        flux_max_bound = config['dfba']['flux_bounds']['max']

        for submodel in model.submodels:
            if are_terms_equivalent(submodel.framework, onto['WC:dynamic_flux_balance_analysis']):
                for rxn in submodel.reactions:
                    if rxn.reversible:
                        flux_min = flux_min_bound_reversible
                    else:
                        flux_min = flux_min_bound_irreversible
                    flux_max = flux_max_bound

                    if rxn.flux_bounds is None:
                        rxn.flux_bounds = wc_lang.core.FluxBounds()

                    if rxn.flux_bounds.min is None or isnan(rxn.flux_bounds.min):
                        rxn.flux_bounds.min = flux_min
                    else:
                        rxn.flux_bounds.min = max(rxn.flux_bounds.min, flux_min)

                    if rxn.flux_bounds.max is None or isnan(rxn.flux_bounds.max):
                        rxn.flux_bounds.max = flux_max
                    else:
                        rxn.flux_bounds.max = min(rxn.flux_bounds.max, flux_max)

                    rxn.flux_bounds.units = unit_registry.parse_units('M s^-1')
        return model
