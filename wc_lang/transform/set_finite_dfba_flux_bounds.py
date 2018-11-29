""" Clip the flux bounds for the reactions in dFBA submodels

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from .core import Transform
from math import isnan
from wc_lang.core import SubmodelAlgorithm
import wc_lang.config.core
config = wc_lang.config.core.get_config()['wc_lang']


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
