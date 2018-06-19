""" Tests of model transforms.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_lang import (Model, Submodel, Reaction, Parameter, SpeciesType, SpeciesTypeType,
                     Species, Compartment, SpeciesCoefficient, RateLawDirection, RateLawEquation, SubmodelAlgorithm)
from wc_lang.transform import MergeAlgorithmicallyLikeSubmodelsTransform
import unittest


class MergeAlgorithmicallyLikeSubmodelsTransformTestCase(unittest.TestCase):
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

        submdl_0 = Submodel(id='submdl_0', algorithm=SubmodelAlgorithm.SSA)
        submdl_1 = Submodel(id='submdl_1', algorithm=SubmodelAlgorithm.SSA)
        submdl_2 = Submodel(id='submdl_2', algorithm=SubmodelAlgorithm.dFBA)
        mdl.submodels.add(submdl_0)
        mdl.submodels.add(submdl_1)
        mdl.submodels.add(submdl_2)

        rxn_0_0 = Reaction(id='rxn_0_0')
        rxn_0_0.participants.add(SpeciesCoefficient(species=specs[0], coefficient=-1))
        rxn_0_0.participants.add(SpeciesCoefficient(species=specs[1], coefficient=-1))
        rxn_0_0.participants.add(SpeciesCoefficient(species=specs[2], coefficient=1))

        rxn_0_1 = Reaction(id='rxn_0_1')
        rxn_0_1.participants.add(SpeciesCoefficient(species=specs[0], coefficient=-1))
        rxn_0_1.participants.add(SpeciesCoefficient(species=specs[1], coefficient=-1))
        rxn_0_1.participants.add(SpeciesCoefficient(species=specs[2], coefficient=1))

        rxn_1_0 = Reaction(id='rxn_1_0')
        rxn_1_0.participants.add(SpeciesCoefficient(species=specs[0], coefficient=-1))
        rxn_1_0.participants.add(SpeciesCoefficient(species=specs[1], coefficient=-1))
        rxn_1_0.participants.add(SpeciesCoefficient(species=specs[3], coefficient=1))

        rxn_1_1 = Reaction(id='rxn_1_1')
        rxn_1_1.participants.add(SpeciesCoefficient(species=specs[0], coefficient=-1))
        rxn_1_1.participants.add(SpeciesCoefficient(species=specs[1], coefficient=-1))
        rxn_1_1.participants.add(SpeciesCoefficient(species=specs[3], coefficient=1))

        rxn_2_0 = Reaction(id='rxn_2_0')
        rxn_2_0.participants.add(SpeciesCoefficient(species=specs[0], coefficient=-1))
        rxn_2_0.participants.add(SpeciesCoefficient(species=specs[1], coefficient=-1))
        rxn_2_0.participants.add(SpeciesCoefficient(species=specs[4], coefficient=1))

        rxn_2_1 = Reaction(id='rxn_2_1')
        rxn_2_1.participants.add(SpeciesCoefficient(species=specs[0], coefficient=-1))
        rxn_2_1.participants.add(SpeciesCoefficient(species=specs[1], coefficient=-1))
        rxn_2_1.participants.add(SpeciesCoefficient(species=specs[4], coefficient=1))

        submdl_0.reactions.add(rxn_0_0)
        submdl_0.reactions.add(rxn_0_1)
        submdl_1.reactions.add(rxn_1_0)
        submdl_1.reactions.add(rxn_1_1)
        submdl_2.reactions.add(rxn_2_0)
        submdl_2.reactions.add(rxn_2_1)

        submdl_0.parameters.create(id='param_0')
        submdl_1.parameters.create(id='param_1')
        submdl_2.parameters.create(id='param_2')

        submdl_0.references.create(id='ref_0')
        submdl_1.references.create(id='ref_1')
        submdl_2.references.create(id='ref_2')

        submdl_0.database_references.create(id='xref_0')
        submdl_1.database_references.create(id='xref_1')
        submdl_2.database_references.create(id='xref_2')

        """ Merge algorithmically-like submodels """
        merged_mdl = mdl.copy()
        MergeAlgorithmicallyLikeSubmodelsTransform().run(merged_mdl)
        merged_submdl_ssa = merged_mdl.submodels.get_one(algorithm=SubmodelAlgorithm.SSA)
        merged_submdl_fba = merged_mdl.submodels.get_one(algorithm=SubmodelAlgorithm.dFBA)

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

        self.assertEqual(len(submdl_0.parameters) + len(submdl_1.parameters), len(merged_submdl_ssa.parameters))
        self.assertEqual(len(submdl_2.parameters), len(merged_submdl_fba.parameters))

        self.assertEqual(len(submdl_0.references) + len(submdl_1.references), len(merged_submdl_ssa.references))
        self.assertEqual(len(submdl_2.references), len(merged_submdl_fba.references))

        self.assertEqual(len(submdl_0.database_references) + len(submdl_1.database_references), len(merged_submdl_ssa.database_references))
        self.assertEqual(len(submdl_2.database_references), len(merged_submdl_fba.database_references))
