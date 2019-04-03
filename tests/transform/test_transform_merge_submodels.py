""" Test that algorithmically-like submodels are correctly merged

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_lang import (Model,
                     Species, SpeciesCoefficient, Reaction,
                     DfbaObjective, DfbaObjectiveExpression, DfbaObjReaction,
                     Evidence)
from wc_lang.transform import MergeAlgorithmicallyLikeSubmodelsTransform
from wc_onto import onto
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
                                          type=onto['WC:metabolite'])
            s = mdl.species.create(species_type=st,
                                   compartment=cmp)
            s.id = s.gen_id()
            species.append(s)

        submdl_0 = mdl.submodels.create(id='submdl_0', framework=onto['WC:stochastic_simulation_algorithm'])
        submdl_1 = mdl.submodels.create(id='submdl_1', framework=onto['WC:stochastic_simulation_algorithm'])
        submdl_2 = mdl.submodels.create(id='submdl_2', framework=onto['WC:dynamic_flux_balance_analysis'])
        submdl_3 = mdl.submodels.create(id='submdl_3', framework=onto['WC:dynamic_flux_balance_analysis'])

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

        mdl.observations.create(id='obs_0')
        mdl.observations.create(id='obs_1')
        mdl.observations.create(id='obs_2')
        mdl.observations.create(id='obs_3')
        mdl.observations.create(id='obs_4')
        mdl.observations.create(id='obs_5')
        mdl.observations.create(id='obs_6')
        mdl.observations.create(id='obs_7')

        ev_0 = Evidence(observation=mdl.observations[0], type=onto['WC:supporting_evidence'])
        ev_1 = Evidence(observation=mdl.observations[1], type=onto['WC:supporting_evidence'])
        ev_2 = Evidence(observation=mdl.observations[2], type=onto['WC:supporting_evidence'])
        ev_3 = Evidence(observation=mdl.observations[3], type=onto['WC:supporting_evidence'])
        submdl_0.evidence.extend([ev_0, ev_1])
        submdl_1.evidence.extend([ev_0, ev_2, ev_3])
        submdl_2.evidence.extend([ev_1, ev_2])
        submdl_3.evidence.extend([ev_3])
        dfba_obj_2.evidence.append(ev_2)
        dfba_obj_3.evidence.append(ev_3)

        mdl.conclusions.create(id='conclusion_0', submodels=[submdl_0, submdl_1],
          evidence=[Evidence(observation=mdl.observations[4], type=onto['WC:supporting_evidence'])])
        mdl.conclusions.create(id='conclusion_1', submodels=[submdl_0, submdl_2],
          evidence=[Evidence(observation=mdl.observations[5], type=onto['WC:supporting_evidence'])])
        mdl.conclusions.create(id='conclusion_2', submodels=[submdl_1, submdl_2], dfba_objs=[dfba_obj_2],
          evidence=[Evidence(observation=mdl.observations[6], type=onto['WC:supporting_evidence'])])
        mdl.conclusions.create(id='conclusion_2', submodels=[submdl_1, submdl_3], dfba_objs=[dfba_obj_3],
          evidence=[Evidence(observation=mdl.observations[7], type=onto['WC:supporting_evidence'])])

        mdl.references.create(id='ref_0', submodels=[submdl_0])
        mdl.references.create(id='ref_1', submodels=[submdl_1])
        mdl.references.create(id='ref_2', submodels=[submdl_2], dfba_objs=[dfba_obj_2])
        mdl.references.create(id='ref_2', submodels=[submdl_3], dfba_objs=[dfba_obj_3])

        submdl_0.identifiers.create(id='xref_0')
        submdl_1.identifiers.create(id='xref_1')
        submdl_2.identifiers.create(id='xref_2', dfba_objs=[dfba_obj_2])
        submdl_3.identifiers.create(id='xref_3', dfba_objs=[dfba_obj_3])

        """ Merge algorithmically-like submodels """
        merged_mdl = mdl.copy()
        MergeAlgorithmicallyLikeSubmodelsTransform().run(merged_mdl)
        merged_submdl_ssa = merged_mdl.submodels.get_one(framework=onto['WC:stochastic_simulation_algorithm'])
        merged_submdl_fba = merged_mdl.submodels.get_one(framework=onto['WC:dynamic_flux_balance_analysis'])

        """ Test submodels merged corrected """
        self.assertEqual(len(merged_mdl.compartments), len(mdl.compartments))
        self.assertEqual(len(merged_mdl.species_types), len(mdl.species_types))
        self.assertEqual(len(merged_mdl.submodels), 2)
        self.assertEqual(len(merged_mdl.parameters), len(mdl.parameters))
        self.assertEqual(len(merged_mdl.reactions), len(mdl.reactions))
        self.assertEqual(len(merged_mdl.dfba_objs), 1)
        self.assertEqual(len(merged_mdl.dfba_obj_reactions), len(mdl.dfba_obj_reactions))
        self.assertEqual(len(merged_mdl.observations), len(mdl.observations))
        self.assertEqual(len(merged_mdl.conclusions), len(mdl.conclusions))
        self.assertEqual(len(merged_mdl.references), len(mdl.references))

        self.assertIn(merged_submdl_ssa.id, [
            '{0}_{1}'.format(submdl_0.id, submdl_1.id),
            '{1}_{0}'.format(submdl_0.id, submdl_1.id),
        ])
        self.assertIn(merged_submdl_fba.id, [
            '{0}_{1}'.format(submdl_2.id, submdl_3.id),
            '{1}_{0}'.format(submdl_2.id, submdl_3.id),
        ])

        self.assertEqual(merged_submdl_ssa.framework, onto['WC:stochastic_simulation_algorithm'])
        self.assertEqual(merged_submdl_fba.framework, onto['WC:dynamic_flux_balance_analysis'])

        self.assertEqual(len(merged_submdl_ssa.get_children(kind='submodel', __type=Species)),
                         len(set(submdl_0.get_children(kind='submodel', __type=Species) +
                                 submdl_1.get_children(kind='submodel', __type=Species))))
        self.assertEqual(len(merged_submdl_fba.get_children(kind='submodel', __type=Species)),
                         len(set(submdl_2.get_children(kind='submodel', __type=Species) +
                                 submdl_3.get_children(kind='submodel', __type=Species))))

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

        self.assertEqual(len(set(submdl_0.conclusions) | set(submdl_1.conclusions)), len(merged_submdl_ssa.conclusions))
        self.assertEqual(len(merged_submdl_fba.conclusions), len(set(submdl_2.conclusions) | set(submdl_3.conclusions)))

        self.assertEqual(len(set(submdl_0.references) | set(submdl_1.references)), len(merged_submdl_ssa.references))
        self.assertEqual(len(merged_submdl_fba.references), len(set(submdl_2.references) | set(submdl_3.references)))

        self.assertEqual(len(set(submdl_0.identifiers) | set(submdl_1.identifiers)), len(merged_submdl_ssa.identifiers))
        self.assertEqual(len(merged_submdl_fba.identifiers), len(set(submdl_2.identifiers) | set(submdl_3.identifiers)))
