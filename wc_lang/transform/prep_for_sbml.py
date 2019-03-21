""" Transform model into SBML-compatible representation and discard information that cannot be exported to SBML.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-03-20
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .core import Transform
from .create_implicit_dist_zero_init_concs import CreateImplicitDistributionZeroInitConcentrationsTransform
from .create_implicit_dfba_ex_rxns import CreateImplicitDfbaExchangeReactionsTransform
from .set_finite_dfba_flux_bounds import SetFiniteDfbaFluxBoundsTransform
from .split_reversible_reactions import SplitReversibleReactionsTransform
from wc_lang.core import Model
import obj_model


class PrepForSbmlTransform(Transform):
    """ Transform model into SBML-compatible representation and discard information that cannot be exported to SBML """

    class Meta(object):
        id = 'PrepForSbml'
        label = 'Transform model into SBML-compatible representation and discard information that cannot be exported to SBML'

    def run(self, model):
        """ Discard information that cannot be exported to SBML

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but only with information that can be exported to SBML
        """
        # transform model into SBML-compatible representation
        transforms = [
            CreateImplicitDistributionZeroInitConcentrationsTransform,
            CreateImplicitDfbaExchangeReactionsTransform,
            SplitReversibleReactionsTransform,
        ]
        for transform in transforms:
            transform().run(model)

        # discard information that can't be exported to SBML
        for obj in model.get_related():
            Meta = obj.__class__.Meta

            assert not set(Meta.sbml_attrs).difference(set(Meta.attributes.keys()))

            discard_attr_names = set(Meta.attributes.keys()).difference(set(Meta.sbml_attrs))
            for attr_name in discard_attr_names:
                attr = Meta.attributes[attr_name]
                none_value = attr.get_none_value()
                setattr(obj, attr_name, none_value)

        return model
