""" Tests of model transforms.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-07-25
:Copyright: 2016, Karr Lab
:License: MIT
"""

from itertools import chain
from wc_lang.core import (Model, Submodel, Reaction, SpeciesType, SpeciesTypeType,
                          Species, Compartment, ReactionParticipant, RateLawDirection, RateLawEquation)
from wc_lang.transform import MergeAlgorithmicallyLikeSubmodelsTransform, SplitReversibleReactionsTransform
from wc_utils.schema.core import RelatedAttribute
import unittest


class TestTransform(unittest.TestCase):
    """ Test model transforms """

    def test_merge_algorithmically_like_submodels(self):
        """ Test that algorithmically-like submodels are correctly merged """

        """ Construct model with 3 submodels: two SSA and one FBA """
        mdl = Model()

        cmp = Compartment(id='comp_0', name='compartment 0')
        mdl.compartments.add(cmp)

        specs = []
        for i in range(5):
            spec_type = SpeciesType(id='spec_type_{}'.format(i), type=SpeciesTypeType.metabolite)
            mdl.species_types.add(spec_type)

            spec = Species(species_type=spec_type, compartment=cmp)
            specs.append(spec)

        submdl_0 = Submodel(id='submdl_0', algorithm='SSA')
        submdl_1 = Submodel(id='submdl_1', algorithm='SSA')
        submdl_2 = Submodel(id='submdl_2', algorithm='FBA')
        mdl.submodels.add(submdl_0)
        mdl.submodels.add(submdl_1)
        mdl.submodels.add(submdl_2)

        rxn_0 = Reaction(id='rxn_0')
        rxn_0.participants.add(ReactionParticipant(species=specs[0], coefficient=-1))
        rxn_0.participants.add(ReactionParticipant(species=specs[1], coefficient=-1))
        rxn_0.participants.add(ReactionParticipant(species=specs[2], coefficient=1))

        rxn_1 = Reaction(id='rxn_1')
        rxn_1.participants.add(ReactionParticipant(species=specs[0], coefficient=-1))
        rxn_1.participants.add(ReactionParticipant(species=specs[1], coefficient=-1))
        rxn_1.participants.add(ReactionParticipant(species=specs[3], coefficient=1))

        rxn_2 = Reaction(id='rxn_2')
        rxn_2.participants.add(ReactionParticipant(species=specs[0], coefficient=-1))
        rxn_2.participants.add(ReactionParticipant(species=specs[1], coefficient=-1))
        rxn_2.participants.add(ReactionParticipant(species=specs[4], coefficient=1))

        submdl_0.reactions.add(rxn_0)
        submdl_1.reactions.add(rxn_1)
        submdl_2.reactions.add(rxn_2)

        """ Merge algorithmically-like submodels """
        merged_mdl = mdl.copy()
        MergeAlgorithmicallyLikeSubmodelsTransform().run(merged_mdl)

        merged_submodels = list(merged_mdl.submodels)

        if merged_submodels[0].algorithm == 'SSA':
            merged_submdl_ssa = merged_submodels[0]
            merged_submdl_fba = merged_submodels[1]
        else:
            merged_submdl_ssa = merged_submodels[1]
            merged_submdl_fba = merged_submodels[0]

        """ Test submodels merged corrected """
        self.assertEqual(len(mdl.compartments), len(merged_mdl.compartments))
        self.assertEqual(len(mdl.species_types), len(merged_mdl.species_types))
        self.assertEqual(2, len(merged_mdl.submodels))
        self.assertEqual(len(mdl.parameters), len(merged_mdl.parameters))
        self.assertEqual(len(mdl.get_reactions()), len(merged_mdl.get_reactions()))

        self.assertIn(merged_submdl_ssa.id, [
            '{0}_{1}'.format(submdl_0.id, submdl_1.id),
            '{1}_{0}'.format(submdl_0.id, submdl_1.id),
        ])
        self.assertEqual(submdl_2.id, merged_submdl_fba.id)

        self.assertEqual(submdl_0.algorithm, merged_submdl_ssa.algorithm)
        self.assertEqual(submdl_2.algorithm, merged_submdl_fba.algorithm)

        self.assertEqual(4, len(merged_submdl_ssa.get_species()))
        self.assertEqual(len(submdl_2.get_species()), len(merged_submdl_fba.get_species()))

        self.assertEqual(len(submdl_0.reactions) + len(submdl_1.reactions), len(merged_submdl_ssa.reactions))
        self.assertEqual(len(submdl_2.reactions), len(merged_submdl_fba.reactions))

    def test_split_reversible_reactions(self):
        model = Model()

        c = model.compartments.create(id='c')

        t0 = model.species_types.create(id='s0', type=SpeciesTypeType.metabolite)
        t1 = model.species_types.create(id='s1', type=SpeciesTypeType.metabolite)
        t2 = model.species_types.create(id='s2', type=SpeciesTypeType.metabolite)

        s0 = Species(species_type=t0, compartment=c)
        s1 = Species(species_type=t1, compartment=c)
        s2 = Species(species_type=t2, compartment=c)

        submodel = model.submodels.create(id='submodel', algorithm='SSA')

        r0 = submodel.reactions.create(id='r0', reversible=True)
        r0.participants.create(species=s0, coefficient=-2)
        r0.participants.create(species=s1, coefficient=3)
        r0_f = r0.rate_laws.create(direction=RateLawDirection.forward, equation=RateLawEquation(expression='a'))
        r0_b = r0.rate_laws.create(direction=RateLawDirection.backward, equation=RateLawEquation(expression='b'))
        r0.references.create(id='ref_0')
        r0.cross_references.create(database='x', id='y')

        r1 = submodel.reactions.create(id='r1', reversible=False)
        r1.participants.create(species=s1, coefficient=-3)
        r1.participants.create(species=s2, coefficient=4)
        r1_f = r1.rate_laws.create(direction=RateLawDirection.forward, equation=RateLawEquation(expression='c'))
        r1_b = r1.rate_laws.create(direction=RateLawDirection.backward, equation=RateLawEquation(expression='d'))
        r1.references.create(id='ref_1')
        r1.cross_references.create(database='xx', id='yy')

        model2 = model.copy()
        submodel2 = model2.submodels.get(id='submodel')
        r0 = submodel2.reactions.get(id='r0')

        SplitReversibleReactionsTransform().run(model2)

        self.assertEqual(set([x.id for x in submodel2.reactions]), set(['r0_forward', 'r0_backward', 'r1']))

        r0_f = submodel2.reactions.get(id='r0_forward')
        r0_b = submodel2.reactions.get(id='r0_backward')
        r1_2 = submodel2.reactions.get(id='r1')
        attr = Reaction.Meta.attributes['participants']
        self.assertEqual(attr.serialize(r0_f.participants), '[c]: (2) s0 ==> (3) s1')
        self.assertEqual(attr.serialize(r0_b.participants), '[c]: (3) s1 ==> (2) s0')
        self.assertEqual(attr.serialize(r1_2.participants), '[c]: (3) s1 ==> (4) s2')
        self.assertEqual(len(r0_f.rate_laws), 1)
        self.assertEqual(len(r0_b.rate_laws), 1)
        self.assertEqual(list(r0_f.rate_laws)[0].direction, RateLawDirection.forward)
        self.assertEqual(list(r0_b.rate_laws)[0].direction, RateLawDirection.forward)
        self.assertEqual(list(r0_f.rate_laws)[0].equation.expression, 'a')
        self.assertEqual(list(r0_b.rate_laws)[0].equation.expression, 'b')
        self.assertEqual(r1_2.rate_laws.get(direction=RateLawDirection.forward).equation.expression, 'c')
        self.assertEqual(r1_2.rate_laws.get(direction=RateLawDirection.backward).equation.expression, 'd')

        self.assertEqual(set([x.id for x in r0_f.references]), set(['ref_0']))
        self.assertEqual(set([x.id for x in r0_b.references]), set(['ref_0']))
        self.assertEqual(set([x.id for x in r1_2.references]), set(['ref_1']))

        self.assertEqual(set([x.id for x in r0_f.cross_references]), set(['y']))
        self.assertEqual(set([x.id for x in r0_b.cross_references]), set(['y']))
        self.assertEqual(set([x.id for x in r1_2.cross_references]), set(['yy']))

        self.assertEqual(r0.submodel, None)
        self.assertEqual(r0.participants, set())
        self.assertEqual(r0.rate_laws, set())
        self.assertEqual(r0.references, set())
        self.assertEqual(r0.cross_references, set())

        for attr_name, attr in chain(Reaction.Meta.attributes.items(), Reaction.Meta.related_attributes.items()):
            if isinstance(attr, RelatedAttribute):
                val = getattr(r0, attr_name)
                self.assertTrue(val is None or (isinstance(val, set) and len(val) == 0))
