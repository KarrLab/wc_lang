""" Test splitting of reversible reactions.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from itertools import chain
from obj_tables import RelatedAttribute
from wc_lang import (Model, Submodel, Reaction, Parameter, SpeciesType, FluxBounds,
                     Species, Compartment, SpeciesCoefficient, RateLawDirection, RateLawExpression,
                     DfbaObjReaction, DfbaObjective, DfbaObjectiveExpression)
from wc_lang.transform import SplitReversibleReactionsTransform
from wc_onto import onto
from wc_utils.util.units import unit_registry
import unittest
import math


class SplitReversibleReactionsTransformTestCase(unittest.TestCase):

    def setUp(self):
        self.model = model = Model()

        c = model.compartments.create(id='c')

        t0 = model.species_types.create(id='s0', type=onto['WC:metabolite'])
        t1 = model.species_types.create(id='s1', type=onto['WC:metabolite'])
        t2 = model.species_types.create(id='s2', type=onto['WC:metabolite'])

        s0 = model.species.create(id='s0[c]', species_type=t0, compartment=c)
        s1 = model.species.create(id='s1[c]', species_type=t1, compartment=c)
        s2 = model.species.create(id='s2[c]', species_type=t2, compartment=c)

        self.submodel = submodel = model.submodels.create(id='submodel', framework='SSA')

        self.r0 = r0 = model.reactions.create(id='r0', reversible=True, submodel=submodel)
        r0.participants.create(species=s0, coefficient=-2)
        r0.participants.create(species=s1, coefficient=3)

        r0_f = r0.rate_laws.create(direction=RateLawDirection.forward, model=model)
        a = model.parameters.create(id='a', value=1., units=unit_registry.parse_units('s^-1'))
        r0_f.expression, error = RateLawExpression.deserialize('a', {Parameter: {'a': a}})
        assert error is None, str(error)

        r0_b = r0.rate_laws.create(direction=RateLawDirection.backward,
                                   expression=RateLawExpression(expression='b'),
                                   model=model)
        b = model.parameters.create(id='b', value=1., units=unit_registry.parse_units('s^-1'))
        r0_b.expression, error = RateLawExpression.deserialize('b', {Parameter: {'b': b}})
        assert error is None, str(error)

        r0.references.create(id='ref_0', model=model)
        r0.identifiers.create(namespace='x', id='y')

        self.r1 = r1 = model.reactions.create(id='r1', reversible=False, submodel=submodel)
        r1.participants.create(species=s1, coefficient=-3)
        r1.participants.create(species=s2, coefficient=4)

        r1_f = r1.rate_laws.create(direction=RateLawDirection.forward, model=model)
        c = model.parameters.create(id='c', value=1., units=unit_registry.parse_units('s^-1'))
        r1_f.expression, error = RateLawExpression.deserialize('c', {Parameter: {'c': c}})
        assert error is None, str(error)

        r1_b = r1.rate_laws.create(direction=RateLawDirection.backward, model=model)
        d = model.parameters.create(id='d', value=1., units=unit_registry.parse_units('s^-1'))
        r1_b.expression, error = RateLawExpression.deserialize('d', {Parameter: {'d': d}})
        assert error is None, str(error)

        r1.references.create(id='ref_1', model=model)
        r1.identifiers.create(namespace='xx', id='yy')

    def test(self):
        submodel2 = self.model.submodels.get_one(id='submodel')
        r0 = submodel2.reactions.get_one(id='r0')

        SplitReversibleReactionsTransform().run(self.model)

        self.assertEqual(set([x.id for x in submodel2.reactions]), set(['r0_forward', 'r0_backward', 'r1']))

        r0_f = submodel2.reactions.get_one(id='r0_forward')
        r0_b = submodel2.reactions.get_one(id='r0_backward')
        r1_2 = submodel2.reactions.get_one(id='r1')
        attr = Reaction.Meta.attributes['participants']
        self.assertEqual(attr.serialize(r0_f.participants), '[c]: (2) s0 ==> (3) s1')
        self.assertEqual(attr.serialize(r0_b.participants), '[c]: (3) s1 ==> (2) s0')
        self.assertEqual(attr.serialize(r1_2.participants), '[c]: (3) s1 ==> (4) s2')
        self.assertEqual(len(r0_f.rate_laws), 1)
        self.assertEqual(len(r0_b.rate_laws), 1)
        self.assertEqual(r0_f.rate_laws[0].direction, RateLawDirection.forward)
        self.assertEqual(r0_b.rate_laws[0].direction, RateLawDirection.forward)
        self.assertEqual(r0_f.rate_laws[0].expression.expression, 'a')
        self.assertEqual(r0_b.rate_laws[0].expression.expression, 'b')
        self.assertEqual(r0_f.rate_laws[0].id, r0_f.rate_laws[0].gen_id())
        self.assertEqual(r0_b.rate_laws[0].id, r0_b.rate_laws[0].gen_id())
        self.assertEqual(r0_f.rate_laws[0].id, 'r0_forward-forward')
        self.assertEqual(r0_b.rate_laws[0].id, 'r0_backward-forward')
        self.assertEqual(r1_2.rate_laws.get_one(direction=RateLawDirection.forward).expression.expression, 'c')
        self.assertEqual(r1_2.rate_laws.get_one(direction=RateLawDirection.backward).expression.expression, 'd')

        self.assertEqual(set([x.id for x in r0_f.references]), set(['ref_0']))
        self.assertEqual(set([x.id for x in r0_b.references]), set(['ref_0']))
        self.assertEqual(set([x.id for x in r1_2.references]), set(['ref_1']))

        self.assertEqual(set([x.id for x in r0_f.identifiers]), set(['y']))
        self.assertEqual(set([x.id for x in r0_b.identifiers]), set(['y']))
        self.assertEqual(set([x.id for x in r1_2.identifiers]), set(['yy']))

        self.assertEqual(r0.submodel, None)
        self.assertEqual(r0.participants, [])
        self.assertEqual(r0.rate_laws, [])
        self.assertEqual(r0.references, [])
        self.assertEqual(r0.identifiers, [])

        for attr_name, attr in chain(Reaction.Meta.attributes.items(), Reaction.Meta.related_attributes.items()):
            if isinstance(attr, RelatedAttribute):
                val = getattr(r0, attr_name)
                self.assertTrue(val is None or (isinstance(val, list) and len(val) == 0))

    def test_dfba_submodel(self):
        self.submodel.framework = onto['WC:dynamic_flux_balance_analysis']

        self.r1.reversible = True
        model = self.model.copy()
        SplitReversibleReactionsTransform().run(model)
        for id in ['r0_forward', 'r0_backward', 'r1_forward', 'r1_backward']:
            rxn = model.reactions.get_one(id=id)
            self.assertEqual(rxn.flux_bounds, None)

        self.r0.flux_bounds = FluxBounds(min=3.,
                                         max=2.,
                                         units=unit_registry.parse_units('M s^-1'))
        model = self.model.copy()
        with self.assertRaisesRegexp(AssertionError, 'min flux bound greater than max'):
            SplitReversibleReactionsTransform().run(model)

        # test NaN bounds
        self.r0.flux_bounds = FluxBounds(units=unit_registry.parse_units('M s^-1'))
        model = self.model.copy()
        SplitReversibleReactionsTransform().run(model)
        r0_f = model.reactions.get_one(id='r0_forward')
        r0_b = model.reactions.get_one(id='r0_backward')
        self.assertEqual(r0_f.flux_bounds.min, 0.)
        self.assertTrue(math.isnan(r0_f.flux_bounds.max))
        self.assertEqual(r0_b.flux_bounds.min, 0.)
        self.assertTrue(math.isnan(r0_b.flux_bounds.max))

        # test typical bounds
        lb, ub = -8, 100
        self.r0.flux_bounds = FluxBounds(min=lb,
                                         max=ub,
                                         units=unit_registry.parse_units('M s^-1'))
        model = self.model.copy()
        SplitReversibleReactionsTransform().run(model)
        r0_f = model.reactions.get_one(id='r0_forward')
        r0_b = model.reactions.get_one(id='r0_backward')
        self.assertEqual(r0_f.flux_bounds.min, 0.)
        self.assertEqual(r0_f.flux_bounds.max, ub)
        self.assertEqual(r0_b.flux_bounds.min, 0.)
        self.assertEqual(r0_b.flux_bounds.max, -lb)

        # test 0 < min bound
        lb, ub = 5, 10
        self.r0.flux_bounds = FluxBounds(min=lb,
                                         max=ub,
                                         units=unit_registry.parse_units('M s^-1'))
        model = self.model.copy()
        SplitReversibleReactionsTransform().run(model)
        r0_f = model.reactions.get_one(id='r0_forward')
        r0_b = model.reactions.get_one(id='r0_backward')
        self.assertEqual(r0_f.flux_bounds.min, lb)
        self.assertEqual(r0_f.flux_bounds.max, ub)
        self.assertEqual(r0_b.flux_bounds.min, 0.)
        self.assertEqual(r0_b.flux_bounds.max, 0.)

        # test max bound < 0
        lb, ub = -50, -10
        self.r0.flux_bounds = FluxBounds(min=lb,
                                         max=ub,
                                         units=unit_registry.parse_units('M s^-1'))
        model = self.model.copy()
        SplitReversibleReactionsTransform().run(model)
        r0_f = model.reactions.get_one(id='r0_forward')
        r0_b = model.reactions.get_one(id='r0_backward')
        self.assertEqual(r0_f.flux_bounds.min, 0.)
        self.assertEqual(r0_f.flux_bounds.max, 0.)
        self.assertEqual(r0_b.flux_bounds.min, -ub)
        self.assertEqual(r0_b.flux_bounds.max, -lb)

    def do_test_dfba_obj_expr(self, obj_expression, reversible_rxns, expected_obj_expr,
                              dfba_obj_reaction_id=None):
        self.submodel.framework = onto['WC:dynamic_flux_balance_analysis']
        model = self.model.copy()
        submodel = model.submodels.get_one(id='submodel')
        reactions = {'r0': model.reactions.get_one(id='r0'),
                     'r1': model.reactions.get_one(id='r1')}
        for rxn in reversible_rxns:
            reactions[rxn].reversible = True
        dfba_rxn = model.dfba_obj_reactions.create(id=dfba_obj_reaction_id, submodel=submodel)
        dfba_obj_expression, error = DfbaObjectiveExpression.deserialize(
            obj_expression, {Reaction: reactions,
                             DfbaObjReaction: {dfba_obj_reaction_id: dfba_rxn}})
        assert error is None, str(error)
        submodel.dfba_obj = DfbaObjective(id='dfba-obj-submodel', expression=dfba_obj_expression)
        error = submodel.dfba_obj.validate()
        assert error is None, str(error)
        SplitReversibleReactionsTransform().run(model)
        self.assertEqual(expected_obj_expr, submodel.dfba_obj.expression.expression)

    def test_dfba_obj_expr(self):
        # test transform dFBA objective expression

        # multiple reactions
        expected_obj_expr = '(r0_forward - r0_backward) + 3. * r1 - (r0_forward - r0_backward) / 2 + 0'
        self.do_test_dfba_obj_expr('r0+ 3.*r1 - r0/2 + 0',
                                   reversible_rxns=['r0'],
                                   expected_obj_expr=expected_obj_expr)

        # multiple reversible reactions
        expected_obj_expr = ("(r0_forward - r0_backward) + 3. * ( r1_forward - r1_backward ) - "
                             "(r0_forward - r0_backward) / 2 + 0")
        self.do_test_dfba_obj_expr('r0+ 3.*r1 - r0/2 + 0',
                                   reversible_rxns=['r0', 'r1'],
                                   expected_obj_expr=expected_obj_expr)

        # a reversible reaction and a DfbaObjReaction with different ids
        expected_obj_expr = ("(r0_forward - r0_backward) + 3. * r1 - "
                             "(r0_forward - r0_backward) / 2 + dfba_rxn")
        self.do_test_dfba_obj_expr('r0+ 3.*r1 - r0/2 + dfba_rxn',
                                   reversible_rxns=['r0'],
                                   expected_obj_expr=expected_obj_expr,
                                   dfba_obj_reaction_id='dfba_rxn')

        # a reversible reaction and a DfbaObjReaction with the same ids
        expected_obj_expr = ("(r0_forward - r0_backward) + 3. * r1 - "
                             "(r0_forward - r0_backward) / 2 + DfbaObjReaction.r0")
        self.do_test_dfba_obj_expr('Reaction.r0+ 3.*r1 - Reaction.r0/2 + DfbaObjReaction.r0',
                                   reversible_rxns=['r0'],
                                   expected_obj_expr=expected_obj_expr,
                                   dfba_obj_reaction_id='r0')
