""" Prepare a model for simulation by wc_sim.

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from .core import Transform
from . import create_implicit_dist_zero_init_concs
from . import create_implicit_dfba_ex_rxns
from . import set_finite_dfba_flux_bounds
from . import split_reversible_reactions


class PrepareForWcSimTransform(Transform):
    """ Prepare a model for simulation by making implicit information in the model
    explicit.

    * Create explicit distributions of initial zero concentrations for species
      that don't have explicit distributions of initial concentrations (i.e.
      implicit distributions of initial zero concentrations)
    * Create implicit exchange reactions for dFBA submodels
    * Clip the flux bounds for the reactions in dFBA submodels to the default
      flux range
    * Split reversible reactions into separate forward and reverse reactions

    Attributes:
        transforms (:obj:`list` of :obj:`Transform`): list of transforms
    """
    DEFAULT_TRANSFORMS = (
        create_implicit_dist_zero_init_concs.CreateImplicitDistributionZeroInitConcentrationsTransform,
        create_implicit_dfba_ex_rxns.CreateImplicitDfbaExchangeReactionsTransform,
        set_finite_dfba_flux_bounds.SetFiniteDfbaFluxBoundsTransform,
        split_reversible_reactions.SplitReversibleReactionsTransform,
    )

    class Meta(object):
        id = 'PrepareForSimulation'
        label = ('Prepare model for simulation by making implicit '
                 'information in the model explicit')

    def __init__(self, transforms=None):
        """
        Args:
            transforms (:obj:`list` of :obj:`Transform`, optional): list of transforms
        """
        self.transforms = transforms or list(self.DEFAULT_TRANSFORMS)

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        for transform in self.transforms:
            transform().run(model)
