""" Prepare a model for simulation.

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from wc_lang.transform.core import Transform
from wc_lang.transform import create_implicit_distribution_zero_init_concentrations
from wc_lang.transform import create_implicit_dfba_ex_rxns
from wc_lang.transform import set_finite_dfba_flux_bounds


class PrepareForWcSimTransform(Transform):
    """ Prepare a model for simulation by making implicit information in the model
    explicit.

    * Create explicit distributions of initial zero concentrations for species
      that don't have explicit distributions of initial concentrations (i.e.
      implicit distributions of initial zero concentrations)
    * Create implicit exchange reactions for dFBA submodels
    * Clip the flux bounds for the reactions in dFBA submodels to the default
      flux range
    """
    TRANSFORMS = (
        create_implicit_distribution_zero_init_concentrations.CreateImplicitDistributionZeroInitConcentrationsTransform,
        create_implicit_dfba_ex_rxns.CreateImplicitDfbaExchangeReactionsTransform,
        set_finite_dfba_flux_bounds.SetFiniteDfbaFluxBoundsTransform,
    )

    class Meta(object):
        id = 'PrepareForSimulation'
        label = ('Prepare model for simulation by making implicit '
                 'information in the model explicit')

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        for transform in self.TRANSFORMS:
            transform().run(model)
