""" Tests of model transforms.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2017-07-25
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang.core import Model, Submodel, Reaction, Species, Compartment, ReactionParticipant
from wc_lang.transform import MergeAlgorithmicallyLikeSubmodels
import unittest


class TestTransform(unittest.TestCase):
    """ Test model transforms """

    def test_merge_algorithmically_like_submodels(self):
        """ Test that algorithmically-like submodels are correctly merged """

        """ Construct model with 3 submodels: two SSA and one FBA """
        cmp = Compartment(id='comp-1', name='compartment 1')

        spec1 = Species(id='spec-1')
        spec2 = Species(id='spec-2')
        spec3 = Species(id='spec-3')
        spec4 = Species(id='spec-4')
        spec5 = Species(id='spec-5')

        submdlA = Submodel(id='submdl-A', algorithm='SSA', species=[spec1, spec2, spec3], reactions=[])
        submdlB = Submodel(id='submdl-B', algorithm='SSA', species=[spec1, spec2, spec4], reactions=[])
        submdlC = Submodel(id='submdl-C', algorithm='FBA', species=[spec1, spec2, spec5], reactions=[])

        rxnA1 = Reaction(id='rxn-A-1', submodel=submdlA, participants=[
            ReactionParticipant(species=spec1, compartment=cmp, coefficient=-1),
            ReactionParticipant(species=spec2, compartment=cmp, coefficient=-1),
            ReactionParticipant(species=spec3, compartment=cmp, coefficient=1),
        ])
        rxnB1 = Reaction(id='rxn-B-1', submodel=submdlB, participants=[
            ReactionParticipant(species=spec1, compartment=cmp, coefficient=-1),
            ReactionParticipant(species=spec2, compartment=cmp, coefficient=-1),
            ReactionParticipant(species=spec4, compartment=cmp, coefficient=1),
        ])
        rxnC1 = Reaction(id='rxn-C-1', submodel=submdlC, participants=[
            ReactionParticipant(species=spec1, compartment=cmp, coefficient=-1),
            ReactionParticipant(species=spec2, compartment=cmp, coefficient=-1),
            ReactionParticipant(species=spec5, compartment=cmp, coefficient=1),
        ])
        submdlA.reactions.append(rxnA1)
        submdlB.reactions.append(rxnB1)
        submdlC.reactions.append(rxnC1)

        mdl = Model(
            compartments=[cmp],
            species=[spec1, spec2, spec3, spec4, spec5],
            submodels=[submdlA, submdlB, submdlC],
            reactions=[rxnA1, rxnB1, rxnC1],
        )

        """ Merge algorithmically-like submodels """
        merged_mdl = MergeAlgorithmicallyLikeSubmodels.transform(mdl)

        """ Test submodels merged corrected """
        self.assertEqual(len(mdl.compartments), len(merged_mdl.compartments))
        self.assertEqual(len(mdl.species), len(merged_mdl.species))
        self.assertEqual(2, len(merged_mdl.submodels))
        self.assertEqual(len(mdl.reactions), len(merged_mdl.reactions))

        self.assertEqual('{0}-{1}'.format(mdl.submodels[0].id, mdl.submodels[1].id), merged_mdl.submodels[0].id)
        self.assertEqual(mdl.submodels[2].id, merged_mdl.submodels[1].id)

        self.assertEqual(mdl.submodels[0].algorithm, merged_mdl.submodels[0].algorithm)
        self.assertEqual(mdl.submodels[2].algorithm, merged_mdl.submodels[1].algorithm)

        self.assertEqual(4, len(merged_mdl.submodels[0].species))
        self.assertEqual(len(mdl.submodels[2].species), len(merged_mdl.submodels[1].species))

        self.assertEqual(
            + len(mdl.submodels[0].reactions) 
            + len(mdl.submodels[1].reactions), len(merged_mdl.submodels[0].reactions))
        self.assertEqual(len(mdl.submodels[2].reactions), len(merged_mdl.submodels[1].reactions))
