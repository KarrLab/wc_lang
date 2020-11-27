""" Transform model into SBML-compatible representation and discard information that cannot be exported to SBML.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2019-03-20
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .core import CompositeTransform
from .create_implicit_dist_zero_init_concs import CreateImplicitDistributionZeroInitConcentrationsTransform
from .create_implicit_dfba_ex_rxns import CreateImplicitDfbaExchangeReactionsTransform
from .nullify_attrs import NullifyAttrsTransform
from .set_finite_dfba_flux_bounds import SetFiniteDfbaFluxBoundsTransform
from .split_reversible_reactions import SplitReversibleReactionsTransform
from wc_lang.core import Model
from wc_onto import onto


class PrepForSbmlTransform(CompositeTransform):
    """ Transform model into SBML-compatible representation and discard information that cannot be exported to SBML """

    dfba_framework = [onto['WC:dynamic_flux_balance_analysis']]
    COMPONENT_TRANSFORMS = (
        CreateImplicitDistributionZeroInitConcentrationsTransform(),
        CreateImplicitDfbaExchangeReactionsTransform(),
        SplitReversibleReactionsTransform(excluded_frameworks=dfba_framework),
        NullifyAttrsTransform('sbml'),
    )

    class Meta(object):
        id = 'PrepForSbml'
        label = 'Transform model into SBML-compatible representation and discard information that cannot be exported to SBML'
