""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang import Model, Reaction, RateLawDirection, FluxBounds
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

    def run(self, model):
        """ Split reversible reactions in submodels into separate forward and backward reactions

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model definition, but with reversible reactions split into separate forward and backward reactions
        """
        for submodel in model.submodels:
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

                    # copy dFBA objective
                    if rxn.dfba_obj_expression:
                        dfba_obj_expr = rxn.dfba_obj_expression
                        parsed_expr = dfba_obj_expr._parsed_expression

                        dfba_obj_expr.expression = parsed_expr.expression = re.sub(
                            r'\b' + rxn.id + r'\b',
                            '({} - {})'.format(rxn_for.id, rxn_bck.id),
                            dfba_obj_expr.expression)

                        parsed_expr._objs[Reaction].pop(rxn.id)
                        parsed_expr._objs[Reaction][rxn_for.id] = rxn_for
                        parsed_expr._objs[Reaction][rxn_bck.id] = rxn_bck
                        parsed_expr.tokenize()

                        rxn.dfba_obj_expression = None
                        rxn_for.dfba_obj_expression = dfba_obj_expr
                        rxn_bck.dfba_obj_expression = dfba_obj_expr

        return model
