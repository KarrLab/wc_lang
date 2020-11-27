""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .core import Transform
from obj_tables.math.expression import ObjTablesTokenCodes
from wc_lang import Model, Reaction, DfbaObjReaction, RateLawDirection, FluxBounds, DfbaObjectiveExpression
from wc_onto import onto
from wc_utils.util.ontology import are_terms_equivalent
import copy
import math
import re


class SplitReversibleReactionsTransform(Transform):
    """ Split reversible reactions in submodels into separate forward and backward reactions """

    class Meta(object):
        id = 'SplitReversibleReactions'
        label = 'Split reversible reactions into separate forward and backward reactions'

    def __init__(self, excluded_frameworks=None):
        self.excluded_frameworks = excluded_frameworks
        """
        Args:
            excluded_frameworks (:obj:`list` of :obj:`pronto.term.Term`, optional): submodels using
                these modeling integration frameworks are not transformed
        """

    def run(self, model):
        """ Split reversible reactions in submodels into separate forward and backward reactions

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but with reversible reactions split into separate
            forward and backward reactions, their flux bounds adjusted accordingly, and the
            dFBA objective expression adjusted accordingly
        """
        for submodel in model.submodels:

            # skip submodels which use an excluded framework
            if self.excluded_frameworks is not None:
                if any([are_terms_equivalent(submodel.framework, excluded_framework)
                        for excluded_framework in self.excluded_frameworks]):
                    continue

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
                        conclusions=rxn.conclusions,
                        identifiers=rxn.identifiers,
                        comments=rxn.comments,
                        references=rxn.references,
                    )
                    rxn_bck = submodel.reactions.create(
                        model=model,
                        id='{}_backward'.format(rxn.id),
                        name='{} (backward)'.format(rxn.name),
                        reversible=False,
                        evidence=rxn.evidence,
                        conclusions=rxn.conclusions,
                        identifiers=rxn.identifiers,
                        comments=rxn.comments,
                        references=rxn.references,
                    )

                    rxn.evidence = []
                    rxn.conclusions = []
                    rxn.identifiers = []
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

                    # copy flux bounds
                    if are_terms_equivalent(submodel.framework, onto['WC:dynamic_flux_balance_analysis']):
                        if rxn.flux_bounds:
                            if not math.isnan(rxn.flux_bounds.min) and not math.isnan(rxn.flux_bounds.max):
                                # assume flux_bounds.min <= flux_bounds.max
                                assert rxn.flux_bounds.min <= rxn.flux_bounds.max, \
                                    f"min flux bound greater than max in {rxn.id}"
                            # Mapping of flux bounds to backward and forward reactions
                            # Principles:
                            # lower bounds must be set, and cannot be negative because that would imply reversible
                            # upper bounds may be NaN if not previously set
                            # the forward reaction uses existing bounds that are positive
                            # the backward reaction uses existing bounds that are negative,
                            # swapping min and max and negating signs

                            # NaNs
                            #                   backward rxn    forward rxn
                            #                   ------------    -----------
                            #   min     max     min     max     min     max
                            #   -----   -----   ----    ----    ----    ----
                            #   NaN     NaN     0       NaN     0       NaN

                            # Values
                            #                   backward rxn    forward rxn
                            #                   ------------    -----------
                            # min/max           min     max     min     max
                            # ---------------   ----    ----    ----    ----
                            # min <= max <= 0  -max     -min    0       0
                            # min <= 0 <= max   0       -min    0       max
                            # 0 <= min <= max   0       0       min     max
                            backward_min = 0.
                            backward_max = float('NaN')
                            forward_min = 0.
                            forward_max = float('NaN')

                            if not math.isnan(rxn.flux_bounds.min):
                                if rxn.flux_bounds.min < 0:
                                    backward_max = -rxn.flux_bounds.min
                                else:
                                    backward_max = 0.
                                    forward_min = rxn.flux_bounds.min

                            if not math.isnan(rxn.flux_bounds.max):
                                if 0 < rxn.flux_bounds.max:
                                    forward_max = rxn.flux_bounds.max
                                else:
                                    forward_max = 0.
                                    backward_min = -rxn.flux_bounds.max

                            rxn_bck.flux_bounds = FluxBounds(min=backward_min,
                                                             max=backward_max,
                                                             units=rxn.flux_bounds.units)
                            rxn_for.flux_bounds = FluxBounds(min=forward_min,
                                                             max=forward_max,
                                                             units=rxn.flux_bounds.units)

                    # transform dFBA objective expression
                    # each dFBA objective expression is transformed for each reaction it uses
                    if rxn.dfba_obj_expression:
                        dfba_obj_expr = rxn.dfba_obj_expression
                        parsed_expr = dfba_obj_expr._parsed_expression

                        # create a new dFBA objective expression
                        # 1. use parsed_expr._obj_tables_tokens to recreate the expression and
                        # the objects it uses, while substituting the split reactions for the reversible reaction
                        new_obj_expr_elements = []
                        all_reactions = {Reaction: {},
                                         DfbaObjReaction: {}}
                        for ot_token in parsed_expr._obj_tables_tokens:
                            if (ot_token.code == ObjTablesTokenCodes.obj_id and
                                issubclass(ot_token.model_type, Reaction) and
                                ot_token.model_id == rxn.id):
                                new_obj_expr_elements.append(f'({rxn_for.id} - {rxn_bck.id})')
                                all_reactions[Reaction][rxn_for.id] = rxn_for
                                all_reactions[Reaction][rxn_bck.id] = rxn_bck
                                continue

                            if (ot_token.code == ObjTablesTokenCodes.obj_id and
                                issubclass(ot_token.model_type, (Reaction, DfbaObjReaction))):
                                new_obj_expr_elements.append(ot_token.token_string)
                                all_reactions[ot_token.model_type][ot_token.model_id] = ot_token.model
                                continue

                            new_obj_expr_elements.append(ot_token.token_string)

                        new_obj_expr = ' '.join(new_obj_expr_elements)

                        # 2. create a new DfbaObjectiveExpression
                        dfba_obj_expr, error = DfbaObjectiveExpression.deserialize(new_obj_expr, all_reactions)
                        assert error is None, str(error)

                        rxn.dfba_obj_expression = None
                        rxn_for.dfba_obj_expression = dfba_obj_expr
                        rxn_bck.dfba_obj_expression = dfba_obj_expr
                        submodel.dfba_obj.expression = dfba_obj_expr

        return model
