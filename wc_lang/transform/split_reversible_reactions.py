""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang import Model, RateLawDirection
import copy

class SplitReversibleReactionsTransform(Transform):
    """ Split reversible reactions into separate forward and backward reactions """

    class Meta(object):
        id = 'SplitReversibleReactions'
        label = 'Split reversible reactions into separate forward and backward reactions'

    def run(self, model):
        """ Split reversible reactions into separate forward and backward reactions

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but with reversible reactions split into separate forward and backward reactions
        """

        for submodel in model.submodels:
            for rxn in submodel.reactions:
                if rxn.reversible:
                    # remove reversible reaction
                    submodel.reactions.remove(rxn)

                    # create separate forward and reverse reactions
                    rxn_for = submodel.reactions.create(
                        id='{}_forward'.format(rxn.id),
                        name='{} (forward)'.format(rxn.name),
                        reversible=False,
                        comments=rxn.comments,
                        references=rxn.references,
                    )
                    rxn_bck = submodel.reactions.create(
                        id='{}_backward'.format(rxn.id),
                        name='{} (backward)'.format(rxn.name),
                        reversible=False,
                        comments=rxn.comments,
                        references=rxn.references,
                    )

                    rxn.references = []

                    # copy participants and negate for backward reaction
                    for part in rxn.participants:
                        rxn_for.participants.append(part)

                        part_back = part.species.species_coefficients.get_one(coefficient=-1 * part.coefficient)
                        if part_back:
                            rxn_bck.participants.append(part_back)
                        else:
                            rxn_bck.participants.create(species=part.species, coefficient=-1 * part.coefficient)

                    rxn.participants = []

                    # copy rate laws
                    law_for = rxn.rate_laws.get_one(direction=RateLawDirection.forward)
                    law_bck = rxn.rate_laws.get_one(direction=RateLawDirection.backward)

                    if law_for:
                        law_for.reaction = rxn_for
                        law_for.direction = RateLawDirection.forward
                        law_for.id = law_for.gen_id(law_for.reaction.id, law_for.direction.name)
                    if law_bck:
                        law_bck.reaction = rxn_bck
                        law_bck.direction = RateLawDirection.forward
                        law_bck.id = law_bck.gen_id(law_bck.reaction.id, law_bck.direction.name)

                    # database references
                    for x_ref in rxn.db_refs:
                        rxn.db_refs.remove(x_ref)
                        rxn_for.db_refs.append(x_ref)
                        rxn_bck.db_refs.append(x_ref)

        return model
