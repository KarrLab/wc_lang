""" Tests of model transforms.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-07-25
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang.core import Model, Submodel, Reaction, SpeciesType, SpeciesTypeType, Species, Compartment, ReactionParticipant
from wc_lang.transform import MergeAlgorithmicallyLikeSubmodelsTransform, SplitReversibleReactionsTransform
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
            spec_type = SpeciesType(id='spec_type_{}'.format(i), type=SpeciesTypeType['metabolite'])
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

    @unittest.skip('me')
    def test_split_reversible_reactions(self):
        pass
