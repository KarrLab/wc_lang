""" Clip the flux bounds for the reactions in dFBA submodels

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from .core import Transform
from math import isnan
from wc_lang.core import SubmodelAlgorithm, ReactionFluxBoundUnit
import wc_lang.config.core


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
        flux_min_bound_reversible = config['dfba']['flux_min_bound_reversible']
        flux_min_bound_irreversible = config['dfba']['flux_min_bound_irreversible']
        flux_max_bound = config['dfba']['flux_max_bound']

        for submodel in model.submodels:
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                for rxn in submodel.reactions:
                    if rxn.reversible:
                        flux_min = flux_min_bound_reversible
                    else:
                        flux_min = flux_min_bound_irreversible
                    flux_max = flux_max_bound

                    if rxn.flux_min is None or isnan(rxn.flux_min):
                        rxn.flux_min = flux_min
                    else:
                        rxn.flux_min = max(rxn.flux_min, flux_min)

                    if rxn.flux_max is None or isnan(rxn.flux_max):
                        rxn.flux_max = flux_max
                    else:
                        rxn.flux_max = min(rxn.flux_max, flux_max)

                    rxn.flux_bound_units = ReactionFluxBoundUnit['M s^-1']
        return model
