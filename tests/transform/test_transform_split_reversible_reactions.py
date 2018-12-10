""" Test splitting of reversible reactions.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from itertools import chain
from wc_lang import (Model, Submodel, Reaction, Parameter, SpeciesType, SpeciesTypeType,
                     Species, Compartment, SpeciesCoefficient, RateLawDirection, RateLawExpression, SubmodelAlgorithm)
from wc_lang.transform import SplitReversibleReactionsTransform
from obj_model import RelatedAttribute
import unittest


class SplitReversibleReactionsTransformTestCase(unittest.TestCase):
    def test(self):
        model = Model()

        c = model.compartments.create(id='c')

        t0 = model.species_types.create(id='s0', type=SpeciesTypeType.metabolite)
        t1 = model.species_types.create(id='s1', type=SpeciesTypeType.metabolite)
        t2 = model.species_types.create(id='s2', type=SpeciesTypeType.metabolite)

        s0 = Species(id='s0[c]', species_type=t0, compartment=c)
        s1 = Species(id='s1[c]', species_type=t1, compartment=c)
        s2 = Species(id='s2[c]', species_type=t2, compartment=c)

        submodel = model.submodels.create(id='submodel', algorithm='SSA')

        r0 = submodel.reactions.create(id='r0', reversible=True)
        r0.participants.create(species=s0, coefficient=-2)
        r0.participants.create(species=s1, coefficient=3)
        r0_f = r0.rate_laws.create(direction=RateLawDirection.forward, expression=RateLawExpression(expression='a'))
        r0_b = r0.rate_laws.create(direction=RateLawDirection.backward, expression=RateLawExpression(expression='b'))
        r0.references.create(id='ref_0')
        r0.db_refs.create(database='x', id='y')

        r1 = submodel.reactions.create(id='r1', reversible=False)
        r1.participants.create(species=s1, coefficient=-3)
        r1.participants.create(species=s2, coefficient=4)
        r1_f = r1.rate_laws.create(direction=RateLawDirection.forward, expression=RateLawExpression(expression='c'))
        r1_b = r1.rate_laws.create(direction=RateLawDirection.backward, expression=RateLawExpression(expression='d'))
        r1.references.create(id='ref_1')
        r1.db_refs.create(database='xx', id='yy')

        model2 = model.copy()
        submodel2 = model2.submodels.get_one(id='submodel')
        r0 = submodel2.reactions.get_one(id='r0')

        SplitReversibleReactionsTransform().run(model2)

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
        self.assertEqual(r0_f.rate_laws[0].id, r0_f.rate_laws[0].gen_id(r0_f.rate_laws[0].reaction.id, r0_f.rate_laws[0].direction.name))
        self.assertEqual(r0_b.rate_laws[0].id, r0_b.rate_laws[0].gen_id(r0_b.rate_laws[0].reaction.id, r0_b.rate_laws[0].direction.name))
        self.assertEqual(r0_f.rate_laws[0].id, 'r0_forward-forward')
        self.assertEqual(r0_b.rate_laws[0].id, 'r0_backward-forward')
        self.assertEqual(r1_2.rate_laws.get_one(direction=RateLawDirection.forward).expression.expression, 'c')
        self.assertEqual(r1_2.rate_laws.get_one(direction=RateLawDirection.backward).expression.expression, 'd')

        self.assertEqual(set([x.id for x in r0_f.references]), set(['ref_0']))
        self.assertEqual(set([x.id for x in r0_b.references]), set(['ref_0']))
        self.assertEqual(set([x.id for x in r1_2.references]), set(['ref_1']))

        self.assertEqual(set([x.id for x in r0_f.db_refs]), set(['y']))
        self.assertEqual(set([x.id for x in r0_b.db_refs]), set(['y']))
        self.assertEqual(set([x.id for x in r1_2.db_refs]), set(['yy']))

        self.assertEqual(r0.submodel, None)
        self.assertEqual(r0.participants, [])
        self.assertEqual(r0.rate_laws, [])
        self.assertEqual(r0.references, [])
        self.assertEqual(r0.db_refs, [])

        for attr_name, attr in chain(Reaction.Meta.attributes.items(), Reaction.Meta.related_attributes.items()):
            if isinstance(attr, RelatedAttribute):
                val = getattr(r0, attr_name)
                self.assertTrue(val is None or (isinstance(val, list) and len(val) == 0))
