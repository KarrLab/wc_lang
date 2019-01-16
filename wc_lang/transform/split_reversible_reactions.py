""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang import Model, Reaction, RateLawDirection
from wc_utils.util.ontology import wcm_ontology
import copy
import pronto
import re


class SplitReversibleReactionsTransform(Transform):
    """ Split reversible reactions in non-dFBA submodels into separate forward and backward reactions """

    class Meta(object):
        id = 'SplitReversibleReactions'
        label = 'Split reversible reactions into separate forward and backward reactions'

    def run(self, model):
        """ Split reversible reactions in non-dFBA submodels into separate forward and backward reactions

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but with reversible reactions split into separate forward and backward reactions
        """
        for submodel in model.submodels:
            if submodel.framework != wcm_ontology['WCM:dynamic_flux_balance_analysis']:
                for rxn in list(submodel.reactions):
                    if rxn.reversible:
                        # remove reversible reaction
                        model.reactions.remove(rxn)
                        submodel.reactions.remove(rxn)

                        # create separate forward and reverse reactions
                        rxn_for = submodel.reactions.create(
                            model=model,
                            id='{}_forward'.format(rxn.id),
                            name='{} (forward)'.format(rxn.name),
                            reversible=False,
                            evidence=rxn.evidence,
                            db_refs=rxn.db_refs,
                            comments=rxn.comments,
                            references=rxn.references,
                        )
                        rxn_bck = submodel.reactions.create(
                            model=model,
                            id='{}_backward'.format(rxn.id),
                            name='{} (backward)'.format(rxn.name),
                            reversible=False,
                            evidence=rxn.evidence,
                            db_refs=rxn.db_refs,
                            comments=rxn.comments,
                            references=rxn.references,
                        )

                        rxn.evidence = []
                        rxn.db_refs = []
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
                            law_for.id = law_for.gen_id()
                        if law_bck:
                            law_bck.reaction = rxn_bck
                            law_bck.direction = RateLawDirection.forward
                            law_bck.id = law_bck.gen_id()

                        # copy dFBA objective: unreachable because only non-dFBA reactions are split
                        if rxn.dfba_obj_expression:
                            dfba_obj_expr = rxn.dfba_obj_expression # pragma: no cover
                            parsed_expr = dfba_obj_expr._parsed_expression # pragma: no cover

                            dfba_obj_expr.expression = parsed_expr.expression = re.sub(
                                r'\b' + rxn.id + r'\b',
                                '({} - {})'.format(rxn_for.id, rxn_bck.id),
                                dfba_obj_expr.expression)  # pragma: no cover

                            parsed_expr._objs[Reaction].pop(rxn.id)  # pragma: no cover
                            parsed_expr._objs[Reaction][rxn_for.id] = rxn_for  # pragma: no cover
                            parsed_expr._objs[Reaction][rxn_bck.id] = rxn_bck  # pragma: no cover
                            parsed_expr.tokenize()  # pragma: no cover

                            rxn.dfba_obj_expression = None  # pragma: no cover
                            rxn_for.dfba_obj_expression = dfba_obj_expr  # pragma: no cover
                            rxn_bck.dfba_obj_expression = dfba_obj_expr  # pragma: no cover

        return model
