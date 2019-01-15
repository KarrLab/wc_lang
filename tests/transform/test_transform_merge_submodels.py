""" Test that algorithmically-like submodels are correctly merged

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_lang import (Model,
                     Species, SpeciesCoefficient, Reaction,
                     DfbaObjective, DfbaObjectiveExpression, DfbaObjReaction)
from wc_lang.transform import MergeAlgorithmicallyLikeSubmodelsTransform
from wc_utils.util.ontology import wcm_ontology
import unittest


class MergeAlgorithmicallyLikeSubmodelsTransformTestCase(unittest.TestCase):
    """ Test that algorithmically-like submodels are correctly merged """

    def test(self):
        """ Test that algorithmically-like submodels are correctly merged """

        """ Construct model with 3 submodels: two SSA and one FBA """
        mdl = Model()

        cmp = mdl.compartments.create(id='comp_0', name='compartment 0')

        species = []
        for i in range(5):
            st = mdl.species_types.create(id='spec_type_{}'.format(i),
                                          type=wcm_ontology['WCM:0000015'])  # metabolite
            s = mdl.species.create(species_type=st,
                                   compartment=cmp)
            s.id = s.gen_id()
            species.append(s)

        submdl_0 = mdl.submodels.create(id='submdl_0', algorithm=wcm_ontology['WCM:0000011'])
        submdl_1 = mdl.submodels.create(id='submdl_1', algorithm=wcm_ontology['WCM:0000011'])
        submdl_2 = mdl.submodels.create(id='submdl_2', algorithm=wcm_ontology['WCM:0000013'])
        submdl_3 = mdl.submodels.create(id='submdl_3', algorithm=wcm_ontology['WCM:0000013'])

        rxn_0_0 = mdl.reactions.create(id='rxn_0_0', submodel=submdl_0)
        rxn_0_0.participants.add(SpeciesCoefficient(species=species[0], coefficient=-1))
        rxn_0_0.participants.add(SpeciesCoefficient(species=species[1], coefficient=-1))
        rxn_0_0.participants.add(SpeciesCoefficient(species=species[2], coefficient=1))

        rxn_0_1 = mdl.reactions.create(id='rxn_0_1', submodel=submdl_0)
        rxn_0_1.participants.add(SpeciesCoefficient(species=species[0], coefficient=-1))
        rxn_0_1.participants.add(SpeciesCoefficient(species=species[1], coefficient=-1))
        rxn_0_1.participants.add(SpeciesCoefficient(species=species[2], coefficient=1))

        rxn_1_0 = mdl.reactions.create(id='rxn_1_0', submodel=submdl_1)
        rxn_1_0.participants.add(SpeciesCoefficient(species=species[0], coefficient=-1))
        rxn_1_0.participants.add(SpeciesCoefficient(species=species[1], coefficient=-1))
        rxn_1_0.participants.add(SpeciesCoefficient(species=species[3], coefficient=1))

        rxn_1_1 = mdl.reactions.create(id='rxn_1_1', submodel=submdl_1)
        rxn_1_1.participants.add(SpeciesCoefficient(species=species[0], coefficient=-1))
        rxn_1_1.participants.add(SpeciesCoefficient(species=species[1], coefficient=-1))
        rxn_1_1.participants.add(SpeciesCoefficient(species=species[3], coefficient=1))

        rxn_2_0 = mdl.reactions.create(id='rxn_2_0', submodel=submdl_2)
        rxn_2_0.participants.add(SpeciesCoefficient(species=species[0], coefficient=-1))
        rxn_2_0.participants.add(SpeciesCoefficient(species=species[1], coefficient=-1))
        rxn_2_0.participants.add(SpeciesCoefficient(species=species[4], coefficient=1))

        rxn_2_1 = mdl.reactions.create(id='rxn_2_1', submodel=submdl_2)
        rxn_2_1.participants.add(SpeciesCoefficient(species=species[0], coefficient=-1))
        rxn_2_1.participants.add(SpeciesCoefficient(species=species[1], coefficient=-1))
        rxn_2_1.participants.add(SpeciesCoefficient(species=species[4], coefficient=1))

        dfba_obj_rxn_2_0 = mdl.dfba_obj_reactions.create(id='dfba_obj_rxn_2_0', submodel=submdl_2)
        dfba_obj_rxn_2_1 = mdl.dfba_obj_reactions.create(id='dfba_obj_rxn_2_1', submodel=submdl_2)
        dfba_obj_rxn_3_0 = mdl.dfba_obj_reactions.create(id='dfba_obj_rxn_3_0', submodel=submdl_3)
        dfba_obj_rxn_3_1 = mdl.dfba_obj_reactions.create(id='dfba_obj_rxn_3_1', submodel=submdl_3)

        dfba_obj_2 = mdl.dfba_objs.create(id='dfba_obj_2', submodel=submdl_2,
                                          expression=DfbaObjectiveExpression.deserialize(
                                              'rxn_2_0 + 2 * rxn_2_1 + 3 * dfba_obj_rxn_2_0 + dfba_obj_rxn_2_1', {
                                                  Reaction: {
                                                      'rxn_2_0': rxn_2_0,
                                                      'rxn_2_1': rxn_2_1,
                                                  },
                                                  DfbaObjReaction: {
                                                      'dfba_obj_rxn_2_0': dfba_obj_rxn_2_0,
                                                      'dfba_obj_rxn_2_1': dfba_obj_rxn_2_1,
                                                  },
                                              })[0])
        dfba_obj_3 = mdl.dfba_objs.create(id='dfba_obj_3', submodel=submdl_3,
                                          expression=DfbaObjectiveExpression.deserialize(
                                              'dfba_obj_rxn_3_0 + dfba_obj_rxn_3_1', {
                                                  Reaction: {
                                                  },
                                                  DfbaObjReaction: {
                                                      'dfba_obj_rxn_3_0': dfba_obj_rxn_3_0,
                                                      'dfba_obj_rxn_3_1': dfba_obj_rxn_3_1,
                                                  }
                                              })[0])

        mdl.parameters.create(id='param_0')
        mdl.parameters.create(id='param_1')
        mdl.parameters.create(id='param_2')

        mdl.evidences.create(id='evidence_0', submodels=[submdl_0, submdl_1])
        mdl.evidences.create(id='evidence_1', submodels=[submdl_0, submdl_2])
        mdl.evidences.create(id='evidence_2', submodels=[submdl_1, submdl_2], dfba_objs=[dfba_obj_2])
        mdl.evidences.create(id='evidence_2', submodels=[submdl_1, submdl_3], dfba_objs=[dfba_obj_3])

        mdl.references.create(id='ref_0', submodels=[submdl_0])
        mdl.references.create(id='ref_1', submodels=[submdl_1])
        mdl.references.create(id='ref_2', submodels=[submdl_2], dfba_objs=[dfba_obj_2])
        mdl.references.create(id='ref_2', submodels=[submdl_3], dfba_objs=[dfba_obj_3])

        submdl_0.db_refs.create(id='xref_0')
        submdl_1.db_refs.create(id='xref_1')
        submdl_2.db_refs.create(id='xref_2', dfba_objs=[dfba_obj_2])
        submdl_3.db_refs.create(id='xref_3', dfba_objs=[dfba_obj_3])

        """ Merge algorithmically-like submodels """
        merged_mdl = mdl.copy()
        MergeAlgorithmicallyLikeSubmodelsTransform().run(merged_mdl)
        merged_submdl_ssa = merged_mdl.submodels.get_one(algorithm=wcm_ontology['WCM:0000011'])
        merged_submdl_fba = merged_mdl.submodels.get_one(algorithm=wcm_ontology['WCM:0000013'])

        """ Test submodels merged corrected """
        self.assertEqual(len(merged_mdl.compartments), len(mdl.compartments))
        self.assertEqual(len(merged_mdl.species_types), len(mdl.species_types))
        self.assertEqual(len(merged_mdl.submodels), 2)
        self.assertEqual(len(merged_mdl.parameters), len(mdl.parameters))
        self.assertEqual(len(merged_mdl.reactions), len(mdl.reactions))
        self.assertEqual(len(merged_mdl.dfba_objs), 1)
        self.assertEqual(len(merged_mdl.dfba_obj_reactions), len(mdl.dfba_obj_reactions))
        self.assertEqual(len(merged_mdl.evidences), len(mdl.evidences))
        self.assertEqual(len(merged_mdl.references), len(mdl.references))

        self.assertIn(merged_submdl_ssa.id, [
            '{0}_{1}'.format(submdl_0.id, submdl_1.id),
            '{1}_{0}'.format(submdl_0.id, submdl_1.id),
        ])
        self.assertIn(merged_submdl_fba.id, [
            '{0}_{1}'.format(submdl_2.id, submdl_3.id),
            '{1}_{0}'.format(submdl_2.id, submdl_3.id),
        ])

        self.assertEqual(merged_submdl_ssa.algorithm, wcm_ontology['WCM:0000011'])
        self.assertEqual(merged_submdl_fba.algorithm, wcm_ontology['WCM:0000013'])

        self.assertEqual(len(merged_submdl_ssa.get_species()), len(set(submdl_0.get_species()) | set(submdl_1.get_species())))
        self.assertEqual(len(merged_submdl_fba.get_species()), len(set(submdl_2.get_species()) | set(submdl_3.get_species())))

        self.assertEqual(len(merged_submdl_ssa.reactions), len(set(submdl_0.reactions) | set(submdl_1.reactions)))
        self.assertEqual(len(merged_submdl_fba.reactions), len(set(submdl_2.reactions) | set(submdl_3.reactions)))

        self.assertEqual(merged_submdl_ssa.dfba_obj, None)
        self.assertEqual(merged_submdl_fba.dfba_obj.id, merged_submdl_fba.dfba_obj.gen_id())

        self.assertEqual(len(merged_submdl_fba.dfba_obj.expression.reactions),
                         len(submdl_2.reactions) + len(submdl_3.reactions))
        self.assertEqual(len(merged_submdl_fba.dfba_obj.expression.dfba_obj_reactions),
                         len(submdl_2.dfba_obj_reactions) + len(submdl_3.dfba_obj_reactions))

        self.assertEqual(merged_submdl_ssa.dfba_obj_reactions, [])
        self.assertEqual(len(merged_submdl_fba.dfba_obj_reactions), len(submdl_2.dfba_obj_reactions) + len(submdl_3.dfba_obj_reactions))

        self.assertEqual(len(set(submdl_0.evidence) | set(submdl_1.evidence)), len(merged_submdl_ssa.evidence))
        self.assertEqual(len(merged_submdl_fba.evidence), len(set(submdl_2.evidence) | set(submdl_3.evidence)))

        self.assertEqual(len(set(submdl_0.references) | set(submdl_1.references)), len(merged_submdl_ssa.references))
        self.assertEqual(len(merged_submdl_fba.references), len(set(submdl_2.references) | set(submdl_3.references)))

        self.assertEqual(len(set(submdl_0.db_refs) | set(submdl_1.db_refs)), len(merged_submdl_ssa.db_refs))
        self.assertEqual(len(merged_submdl_fba.db_refs), len(set(submdl_2.db_refs) | set(submdl_3.db_refs)))
