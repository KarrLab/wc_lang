""" Prepare a model for simulation.

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from wc_lang.transform.core import Transform
from wc_lang.transform.create_implicit_zero_concentrations import CreateImplicitZeroConcentrationsTransform
from wc_lang.transform.create_implicit_dfba_ex_rxns import CreateImplicitDfbaExchangeReactionsTransform
from wc_lang.transform.set_finite_dfba_flux_bounds import SetFiniteDfbaFluxBoundsTransform


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
