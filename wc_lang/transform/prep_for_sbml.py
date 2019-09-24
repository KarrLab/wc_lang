""" Transform model into SBML-compatible representation and discard information that cannot be exported to SBML.

:Author: Jonathan Karr <karr@mssm.edu>
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
import obj_tables


class PrepForSbmlTransform(CompositeTransform):
    """ Transform model into SBML-compatible representation and discard information that cannot be exported to SBML """

    COMPONENT_TRANSFORMS = (
        CreateImplicitDistributionZeroInitConcentrationsTransform(),
        CreateImplicitDfbaExchangeReactionsTransform(),
        SplitReversibleReactionsTransform(),
        NullifyAttrsTransform('sbml'),
    )

    class Meta(object):
        id = 'PrepForSbml'
        label = 'Transform model into SBML-compatible representation and discard information that cannot be exported to SBML'
