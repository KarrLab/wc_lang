""" Tests of core

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2016-11-10
:Copyright: 2016-2018, Karr Lab
:License: MIT
"""

import math
import os
import pytest
import re
import unittest
import wc_lang
import wc_lang.config.core
from obj_tables import InvalidAttribute
from obj_tables.expression import ExpressionManyToOneAttribute
from test.support import EnvironmentVarGuard
from wc_lang.core import (Model, Taxon, TaxonRank, Submodel,
                          DfbaObjective, DfbaObjectiveExpression,
                          Reaction, Compartment, InitVolume,
                          SpeciesType, Species,
                          SpeciesCoefficient, Parameter, Reference,
                          Identifier,
                          RateLaw, RateLawExpression, RateLawDirection,
                          Function, FunctionExpression,
                          Observable, ObservableExpression,
                          StopCondition, StopConditionExpression,
                          DistributionInitConcentration, DfbaObjSpecies, DfbaObjReaction,
                          FluxBounds,
                          Observation, ObservationEnv, Evidence, EvidenceManyToManyAttribute, Conclusion, Author, Change,
                          ReactionParticipantAttribute, Expression,
                          InvalidObject, Validator,
                          ChemicalStructure, ChemicalStructureFormat, ChemicalStructureAlphabet)
from wc_lang.io import Reader
from wc_utils.util.chem import EmpiricalFormula
from wc_onto import onto
from wc_utils.util.units import unit_registry


class TestCore(unittest.TestCase):

    def setUp(self):
        self.model = mdl = Model(id='model', name='test model', version='0.0.1', wc_lang_version='0.0.1')

        mdl.taxon = Taxon(id='taxon', name='test taxon', rank=TaxonRank.species)

        self.comp_0 = comp_0 = mdl.compartments.create(id='comp_0', name='compartment 0')
        self.comp_1 = comp_1 = mdl.compartments.create(id='comp_1', name='compartment 1')
        self.compartments = compartments = [comp_0, comp_1]

        self.species_types = species_types = []
        self.species = species = []
        self.distribution_init_concentrations = distribution_init_concentrations = []
        for i in range(8):
            spec_type = mdl.species_types.create(
                id='spec_type_{}'.format(i),
                name='species type {}'.format(i),
                type=onto['WC:metabolite'],
                structure=ChemicalStructure(
                    value='C' * i + 'H' * (i + 1),
                    empirical_formula=EmpiricalFormula('C{}H{}'.format(i, i + 1)),
                    molecular_weight=12 * (i + 1),
                    charge=i + 1))
            species_types.append(spec_type)

            if i != 3:
                spec = Species(species_type=spec_type, compartment=comp_0)
            else:
                spec = Species(species_type=spec_type, compartment=comp_1)
            spec.id = spec.gen_id()
            spec.model = mdl
            species.append(spec)

            conc = DistributionInitConcentration(species=spec, mean=3 * i)
            conc.id = conc.gen_id()
            conc.model = mdl
            distribution_init_concentrations.append(conc)

        objects = {
            Species: {
                'spec_type_0[comp_0]': species[0],
                'spec_type_1[comp_0]': species[1],
                'spec_type_2[comp_0]': species[2],
                'spec_type_3[comp_1]': species[3],
            },
        }
        obs_0 = mdl.observables.create(id='obs_0')
        obs_1 = mdl.observables.create(id='obs_1')
        obs_0.expression, _ = ObservableExpression.deserialize('spec_type_0[comp_0] + spec_type_1[comp_0]', objects)
        obs_1.expression, _ = ObservableExpression.deserialize('spec_type_1[comp_0] + spec_type_2[comp_1]', objects)

        objects = {
            Species: {
                'spec_type_0[comp_0]': species[0],
                'spec_type_1[comp_0]': species[1],
                'spec_type_2[comp_0]': species[2],
                'spec_type_3[comp_1]': species[3],
            },
        }
        func_0 = mdl.functions.create(id='func_0')
        func_1 = mdl.functions.create(id='func_1')
        func_0.expression, _ = FunctionExpression.deserialize('spec_type_0[comp_0] + spec_type_1[comp_0]', objects)
        func_1.expression, _ = FunctionExpression.deserialize('spec_type_1[comp_0] + spec_type_2[comp_1]', objects)

        self.dfba_obj_reaction = dfba_obj_reaction = DfbaObjReaction(
            id='dfba_obj_reaction_1',
            name='dFBA objective reaction',
            model=mdl,
            comments="Nobody will ever deprive the American people of the right to vote except the "
            "American people themselves")
        DfbaObjReaction.get_manager().insert_all_new()

        dfba_obj_species = []
        for i in range(2):
            tmp = dfba_obj_reaction.dfba_obj_species.create(
                value=2 * (float(i) - 0.5),  # create a reactant and a product
                species=species[i])
            tmp.id = tmp.gen_id()
            dfba_obj_species.append(tmp)
        self.dfba_obj_species = dfba_obj_species

        self.submdl_0 = submdl_0 = mdl.submodels.create(
            id='submodel_0', name='submodel 0', framework=onto['WC:stochastic_simulation_algorithm'])
        self.submdl_1 = submdl_1 = mdl.submodels.create(
            id='submodel_1', name='submodel 1', framework=onto['WC:stochastic_simulation_algorithm'])
        self.submdl_2 = submdl_2 = mdl.submodels.create(
            id='submodel_2', name='submodel 2', framework=onto['WC:dynamic_flux_balance_analysis'],
            dfba_obj_reactions=[dfba_obj_reaction])
        self.submodels = submodels = [submdl_0, submdl_1, submdl_2]

        self.parameters = parameters = []
        for i in range(3):
            param = mdl.parameters.create(
                id='param_{}'.format(i), name='parameter {}'.format(i),
                value=i * 4, units=unit_registry.parse_units('dimensionless'))
            param.submodels = submodels[i:i + 1]
            parameters.append(param)

        self.rxn_0 = rxn_0 = submdl_0.reactions.create(
            id='rxn_0', name='reaction 0', model=mdl)
        rxn_0.participants.create(species=species[0], coefficient=-2)
        rxn_0.participants.create(species=species[1], coefficient=-3.5)
        rxn_0.participants.create(species=species[2], coefficient=1)
        expression = RateLawExpression(
            expression='k_cat_0 * {0} / (k_m_0 + {0})'.format(species[5].serialize()),
            species=species[5:6])
        parameters.append(expression.parameters.create(id='k_cat_0', value=2, units=unit_registry.parse_units('s^-1'),
                                                       model=mdl))
        parameters.append(expression.parameters.create(id='k_m_0', value=1, units=unit_registry.parse_units('molecule'),
                                                       model=mdl))
        rate_law_0 = rxn_0.rate_laws.create(
            model=mdl,
            direction=RateLawDirection.forward,
            expression=expression)
        rate_law_0.id = rate_law_0.gen_id()

        self.rxn_1 = rxn_1 = submdl_1.reactions.create(
            id='rxn_1', name='reaction 1', model=mdl)
        rxn_1.participants.create(species=species[0], coefficient=-2)
        rxn_1.participants.create(species=species[1], coefficient=-3)
        rxn_1.participants.create(species=species[3], coefficient=2)
        expression = RateLawExpression(
            expression='k_cat_1 * {0} / (k_m_1 + {0})'.format(species[6].serialize()),
            species=species[6:7])
        parameters.append(expression.parameters.create(id='k_cat_1', value=2, units=unit_registry.parse_units('s^-1'),
                                                       model=mdl))
        parameters.append(expression.parameters.create(id='k_m_1', value=1, units=unit_registry.parse_units('molecule'),
                                                       model=mdl))
        rate_law_1 = rxn_1.rate_laws.create(
            model=mdl,
            direction=RateLawDirection.forward,
            expression=expression)
        rate_law_1.id = rate_law_1.gen_id()

        self.rxn_2 = rxn_2 = submdl_2.reactions.create(
            id='rxn_2', name='reaction 2', model=mdl)
        rxn_2.participants.create(species=species[0], coefficient=-2)
        rxn_2.participants.create(species=species[1], coefficient=-3)
        rxn_2.participants.create(species=species[4], coefficient=1)
        expression = RateLawExpression(
            expression='{1} * {0} / ({2} + {0})'.format(species[7].serialize(), parameters[0].id, parameters[1].id),
            species=species[7:8],
            parameters=parameters[0:2])
        rate_law_2 = rxn_2.rate_laws.create(
            model=mdl,
            direction=RateLawDirection.forward,
            expression=expression)
        rate_law_2.id = rate_law_2.gen_id()

        Reaction.get_manager().insert_all_new()

        self.reactions = [rxn_0, rxn_1, rxn_2]
        self.rate_laws = [rate_law_0, rate_law_1, rate_law_2]

        self.dfba_obj = of = DfbaObjective()
        of.model = mdl
        of.submodel = submdl_2
        of.expression = DfbaObjectiveExpression()
        of.expression.reactions.append(rxn_1)
        of.expression.reactions.append(rxn_2)
        dfba_obj_reaction.dfba_obj_expression = of.expression

        stop_cond_1 = mdl.stop_conditions.create(id='stop_cond_1')
        stop_cond_1.expression, _ = StopConditionExpression.deserialize('2 > 1', {})

        mdl.observations.create(id='obs_0')
        mdl.observations.create(id='obs_2')

        mdl.observations[0].evidence.create(quality=10.)
        mdl.observations[1].evidence.create(quality=20.)

        mdl.conclusions.create(id='int_0', submodels=[submdl_0])
        mdl.conclusions.create(id='int_2', submodels=[submdl_2])

        self.references = references = []
        self.identifiers = identifiers = []
        for i in range(3):
            ref = parameters[i].references.create(
                id='ref_{}'.format(i), name='reference {}'.format(i),
                model=mdl,
                type=None)
            references.append(ref)

            x_ref = ref.identifiers.create(namespace='x', id='y' * (i + 1))
            identifiers.append(x_ref)

        mdl.authors.create(id='First_Last', name='First Last', last_name='Last', first_name='First')
        mdl.authors.create(id='First_Middle_Last', name='First Middel Last', last_name='Last', first_name='First', middle_name='Middle')
        mdl.changes.create(id='change_0')
        mdl.changes.create(id='change_1')

    def test_default_wc_lang_version(self):
        model = Model()
        self.assertEqual(model.wc_lang_version, wc_lang.__version__)

        model = Model(wc_lang_version='xxx')
        self.assertEqual(model.wc_lang_version, 'xxx')

        model = Model(revision='xxx')
        self.assertEqual(model.revision, 'xxx')

    def test_reverse_references(self):
        mdl = self.model

        self.assertEqual(set(mdl.submodels), set(self.submodels))
        self.assertEqual(set(mdl.compartments), set(self.compartments))
        self.assertEqual(set(mdl.species_types), set(self.species_types))
        self.assertEqual(set(mdl.parameters), set(self.parameters))

        # submodel
        for reaction, submodel in zip(self.reactions, self.submodels):
            self.assertEqual(submodel.reactions, [reaction])

        for submodel in self.submodels:
            self.assertEqual(submodel.identifiers, [])
            self.assertEqual(submodel.references, [])

        # compartment
        self.assertEqual(set(self.compartments[0].species), set(self.species[0:3] + self.species[4:]))
        self.assertEqual(self.compartments[1].species, self.species[3:4])

        for compartment in self.compartments:
            self.assertEqual(compartment.identifiers, [])
            self.assertEqual(compartment.references, [])

        # species type
        for species_type, species in zip(self.species_types, self.species):
            self.assertEqual(species_type.species, [species])

        for species_type in self.species_types:
            self.assertEqual(species_type.identifiers, [])
            self.assertEqual(species_type.references, [])

        # specie
        for species_type, species in zip(self.species_types, self.species):
            self.assertEqual(species.species_type, species_type)

        for i in range(len(self.species)):
            if i != 3:
                self.assertEqual(self.species[i].compartment, self.compartments[0])
            else:
                self.assertEqual(self.species[i].compartment, self.compartments[1])

        self.assertEqual(len(self.species[0].species_coefficients), 3)
        self.assertEqual(len(self.species[1].species_coefficients), 3)
        self.assertEqual(len(self.species[2].species_coefficients), 1)
        self.assertEqual(len(self.species[3].species_coefficients), 1)
        self.assertEqual(len(self.species[4].species_coefficients), 1)
        self.assertEqual(len(self.species[5].species_coefficients), 0)
        self.assertEqual(len(self.species[6].species_coefficients), 0)
        self.assertEqual(len(self.species[7].species_coefficients), 0)

        self.assertEqual(len(self.species[0].rate_law_expressions), 0)
        self.assertEqual(len(self.species[1].rate_law_expressions), 0)
        self.assertEqual(len(self.species[2].rate_law_expressions), 0)
        self.assertEqual(len(self.species[3].rate_law_expressions), 0)
        self.assertEqual(len(self.species[4].rate_law_expressions), 0)
        self.assertEqual(len(self.species[5].rate_law_expressions), 1)
        self.assertEqual(len(self.species[6].rate_law_expressions), 1)
        self.assertEqual(len(self.species[7].rate_law_expressions), 1)

        # reaction
        for reaction, submodel in zip(self.reactions, self.submodels):
            self.assertEqual(reaction.submodel, submodel)

        self.assertEqual(set(x.species for x in self.reactions[0].participants),
                         set([self.species[0], self.species[1], self.species[2]]))
        self.assertEqual(set(x.species for x in self.reactions[1].participants),
                         set([self.species[0], self.species[1], self.species[3]]))
        self.assertEqual(set(x.species for x in self.reactions[2].participants),
                         set([self.species[0], self.species[1], self.species[4]]))

        self.assertEqual(self.reactions[0].rate_laws[0].expression.species, self.species[5:6])
        self.assertEqual(self.reactions[1].rate_laws[0].expression.species, self.species[6:7])
        self.assertEqual(self.reactions[2].rate_laws[0].expression.species, self.species[7:8])

        for reaction in self.reactions:
            self.assertEqual(reaction.identifiers, [])
            self.assertEqual(reaction.references, [])
            self.assertEqual(len(reaction.rate_laws), 1)

        # dFBA objective species
        for i in range(len(self.dfba_obj_species)):
            # submodels
            self.assertEqual(self.dfba_obj_species[i].dfba_obj_reaction, self.dfba_obj_reaction)
            # self.assertEqual(self.dfba_obj_reaction.submodels[0], self.submodels[2])
            # species types
            self.assertEqual(self.dfba_obj_species[i].species, self.species[i])
            self.assertEqual(self.dfba_obj_species[i], self.species[i].dfba_obj_species[0])

        # parameters
        for reference, parameter in zip(self.references, self.parameters):
            self.assertEqual(parameter.references, [reference])

        for parameter in self.parameters:
            self.assertEqual(parameter.model, mdl)

        # references
        for reference, parameter in zip(self.references, self.parameters):
            self.assertEqual(reference.parameters, [parameter])
            self.assertEqual(parameter.references, [reference])

        for reference, identifier in zip(self.references, self.identifiers):
            self.assertEqual(reference.identifiers, [identifier])
            self.assertEqual(identifier.references, [reference])

        # reaction participant
        for species in self.species[0:5]:
            self.assertEqual(set(x.species for x in species.species_coefficients), set([species]))

        for reaction in self.reactions:
            for part in reaction.participants:
                self.assertIn(reaction, part.reactions)
            self.assertEqual(set(x.reaction for x in reaction.rate_laws), set([reaction]))

        # identifiers
        for reference, identifier in zip(self.references, self.identifiers):
            self.assertEqual(reference.identifiers, [identifier])
            self.assertEqual(identifier.references, [reference])

    def test_taxon_rank_class(self):
        self.assertEqual(TaxonRank['class'], TaxonRank['classis'])
        self.assertEqual(TaxonRank.__getattr__('class'), TaxonRank['classis'])

    def test_model_get_species(self):
        model = self.model
        self.assertEqual(set(model.get_species()), set(self.species))
        self.assertNotEqual(set(model.get_species()), set())
        self.assertEqual(set(model.get_species(__type=Species)), set(self.species))
        self.assertEqual(model.get_species(__type=Model), [])

    def test_model_get_distribution_init_concentrations(self):
        model = self.model
        self.assertEqual(set(model.get_distribution_init_concentrations()), set(model.distribution_init_concentrations))
        self.assertNotEqual(set(model.get_distribution_init_concentrations()), set())
        self.assertEqual(set(model.get_distribution_init_concentrations(__type=DistributionInitConcentration)),
                         set(model.distribution_init_concentrations))
        self.assertEqual(model.get_distribution_init_concentrations(__type=Model), [])

    def test_model_get_observables(self):
        model = self.model
        self.assertEqual(set(model.get_observables()), set(model.observables))
        self.assertNotEqual(set(model.get_observables()), set())
        self.assertEqual(set(model.get_observables(__type=Observable)), set(model.observables))
        self.assertEqual(model.get_observables(__type=Model), [])

    def test_model_get_functions(self):
        model = self.model
        self.assertEqual(set(model.get_functions()), set(model.functions))
        self.assertNotEqual(set(model.get_functions()), set())
        self.assertEqual(set(model.get_functions(__type=Function)), set(model.functions))
        self.assertEqual(model.get_functions(__type=Model), [])

    def test_model_get_rate_laws(self):
        model = self.model
        self.assertEqual(set(model.get_rate_laws()), set(model.rate_laws))
        self.assertNotEqual(set(model.get_rate_laws()), set())
        self.assertEqual(set(model.get_rate_laws(__type=RateLaw)), set(model.rate_laws))
        self.assertEqual(model.get_rate_laws(__type=Model), [])

    def test_model_get_dfba_objs(self):
        model = self.model
        self.assertEqual(set(model.get_dfba_objs()), set(model.dfba_objs))
        self.assertNotEqual(set(model.get_dfba_objs()), set())
        self.assertEqual(set(model.get_dfba_objs(__type=DfbaObjective)), set(model.dfba_objs))
        self.assertEqual(model.get_dfba_objs(__type=Model), [])

    def test_model_get_dfba_obj_reactions(self):
        model = self.model
        self.assertEqual(set(model.get_dfba_obj_reactions()), set(model.dfba_obj_reactions))
        self.assertNotEqual(set(model.get_dfba_obj_reactions()), set())
        self.assertEqual(set(model.get_dfba_obj_reactions(__type=DfbaObjReaction)), set(model.dfba_obj_reactions))
        self.assertEqual(model.get_dfba_obj_reactions(__type=Model), [])

    def test_model_get_stop_conditions(self):
        model = self.model
        self.assertEqual(set(model.get_stop_conditions()), set(model.stop_conditions))
        self.assertNotEqual(set(model.get_stop_conditions()), set())
        self.assertEqual(set(model.get_stop_conditions(__type=StopCondition)), set(model.stop_conditions))
        self.assertEqual(model.get_stop_conditions(__type=Model), [])

    def test_model_get_observations(self):
        model = self.model
        self.assertEqual(set(model.get_observations()), set(model.observations))
        self.assertNotEqual(set(model.get_observations()), set())
        self.assertEqual(set(model.get_observations(__type=Observation)), set(model.observations))
        self.assertEqual(model.get_observations(__type=Model), [])

    def test_model_get_evidence(self):
        model = self.model

        evidence = []
        for obs in model.observations:
            evidence.extend(obs.evidence)

        self.assertEqual(set(model.get_evidence()), set(evidence))
        self.assertNotEqual(set(model.get_evidence()), set())
        self.assertEqual(set(model.get_evidence(__type=Evidence)), set(evidence))
        self.assertEqual(model.get_evidence(__type=Model), [])

    def test_model_get_conclusions(self):
        model = self.model
        self.assertEqual(set(model.get_conclusions()), set(model.conclusions))
        self.assertNotEqual(set(model.get_conclusions()), set())
        self.assertEqual(set(model.get_conclusions(__type=Conclusion)), set(model.conclusions))
        self.assertEqual(model.get_conclusions(__type=Model), [])

    def test_model_get_authors(self):
        model = self.model
        self.assertEqual(set(model.get_authors()), set(model.authors))
        self.assertNotEqual(set(model.get_authors()), set())
        self.assertEqual(set(model.get_authors(__type=Author)), set(model.authors))
        self.assertEqual(model.get_authors(__type=Model), [])

    def test_model_get_changes(self):
        model = self.model
        self.assertEqual(set(model.get_changes()), set(model.changes))
        self.assertNotEqual(set(model.get_changes()), set())
        self.assertEqual(set(model.get_changes(__type=Change)), set(model.changes))
        self.assertEqual(model.get_changes(__type=Model), [])

    def test_submodel_validate(self):
        submdl = Submodel(id='submodel')
        self.assertEqual(submdl.validate(), None)

        submdl = Submodel(id='sub-model')
        self.assertNotEqual(submdl.validate(), None)

    def test_submodel_get_compartments(self):
        compartments = self.compartments

        self.assertEqual(self.submdl_0.get_children(kind='submodel', __type=Compartment), compartments[0:1])
        self.assertEqual(set(self.submdl_1.get_children(kind='submodel', __type=Compartment)), set(compartments))
        self.assertEqual(set(self.submdl_2.get_children(kind='submodel', __type=Compartment)), set(compartments))

    def test_submodel_get_species_types(self):
        species_types = self.species_types
        self.assertEqual(set(self.submdl_0.get_children(kind='submodel', __type=SpeciesType)), set([
            species_types[0], species_types[1], species_types[2], species_types[5],
        ]))
        self.assertEqual(set(self.submdl_1.get_children(kind='submodel', __type=SpeciesType)), set([
            species_types[0], species_types[1], species_types[3], species_types[6],
        ]))
        self.assertEqual(set(self.submdl_2.get_children(kind='submodel', __type=SpeciesType)), set([
            species_types[0], species_types[1], species_types[3], species_types[4], species_types[6], species_types[7],
        ]))

    def test_submodel_get_species(self):
        species = self.species
        self.assertEqual(set(self.submdl_0.get_children(kind='submodel', __type=Species)), set([
            species[0], species[1], species[2], species[5],
        ]))
        self.assertEqual(set(self.submdl_1.get_children(kind='submodel', __type=Species)), set([
            species[0], species[1], species[3], species[6],
        ]))
        self.assertEqual(set(self.submdl_2.get_children(kind='submodel', __type=Species)), set([
            species[0], species[1], species[3], species[4], species[6], species[7],
        ]))

    def test_submodel_get_observation(self):
        species = self.species
        conc = species[2].conclusions.create()
        self.assertEqual(set(self.submdl_0.get_children(kind='submodel', __type=Conclusion)),
                         set([conc]) | set(self.model.conclusions[0:1]))
        self.assertEqual(self.submdl_1.get_children(kind='submodel', __type=Conclusion), [])
        self.assertEqual(set(self.submdl_2.get_children(kind='submodel', __type=Conclusion)),
                         set(self.model.conclusions[1:2]))

    def test_submodel_get_conclusions(self):
        species = self.species
        conclusion = species[2].conclusions.create()
        self.assertEqual(set(self.submdl_0.get_children(kind='submodel', __type=Conclusion)),
                         set([conclusion]) | set(self.model.conclusions[0:1]))
        self.assertEqual(self.submdl_1.get_children(kind='submodel', __type=Conclusion), [])
        self.assertEqual(set(self.submdl_2.get_children(kind='submodel', __type=Conclusion)),
                         set(self.model.conclusions[1:2]))

    def test_submodel_get_references(self):
        species = self.species
        species[0].references = self.references[0:1]
        species[2].references = self.references[1:2]
        species[3].references = self.references[2:3]
        self.assertEqual(set(self.submdl_0.get_children(kind='submodel', __type=Reference)),
                         set([self.references[0], self.references[1]]))
        self.assertEqual(set(self.submdl_1.get_children(kind='submodel', __type=Reference)),
                         set([self.references[0], self.references[2]]))
        self.assertEqual(set(self.submdl_2.get_children(kind='submodel', __type=Reference)),
                         set(self.references[0:3]))

    def test_submodel_get_components(self):
        self.assertIn(self.submdl_0.reactions[0], self.submdl_0.get_children(kind='submodel'))
        self.assertIn(self.submdl_1.reactions[0], self.submdl_1.get_children(kind='submodel'))
        self.assertIn(self.submdl_1.reactions[0], self.submdl_2.get_children(kind='submodel'))

    def test_compartment_validate(self):
        comp = Compartment(id='c', geometry=onto['WC:3D_compartment'],
                           init_density=Parameter(units=unit_registry.parse_units('g l^-1')))
        self.assertEqual(comp.validate(), None)

        comp = Compartment(id='', geometry=onto['WC:3D_compartment'],
                           init_density=Parameter(units=unit_registry.parse_units('g l^-1')))
        self.assertNotEqual(comp.validate(), None)

        comp = Compartment(id='c', geometry=onto['WC:3D_compartment'])
        self.assertNotEqual(comp.validate(), None)

        comp = Compartment(id='c', geometry=onto['WC:3D_compartment'], init_density=Parameter(units=None))
        self.assertNotEqual(comp.validate(), None)

    def test_reaction_get_species(self):
        species = self.species
        self.assertEqual(set(self.rxn_0.get_species()), set([
            species[0], species[1], species[2], species[5],
        ]))
        self.assertEqual(set(self.rxn_1.get_species()), set([
            species[0], species[1], species[3], species[6],
        ]))
        self.assertEqual(set(self.rxn_2.get_species()), set([
            species[0], species[1], species[4], species[7],
        ]))

        self.assertEqual(set(self.rxn_0.get_species(__type=Species)), set([
            species[0], species[1], species[2], species[5],
        ]))
        self.assertEqual(self.rxn_0.get_species(__type=Model), [])

    def test_get_components(self):
        mdl = self.model

        self.assertEqual(set(mdl.get_compartments()), set(self.compartments))
        self.assertEqual(set(mdl.get_compartments(__type=Compartment)), set(self.compartments))
        self.assertEqual(set(mdl.get_compartments(__type=Submodel)), set())

        self.assertEqual(set(mdl.get_species_types()), set(self.species_types))
        self.assertEqual(set(mdl.get_species_types(__type=SpeciesType)), set(self.species_types))
        self.assertEqual(set(mdl.get_species_types(__type=Compartment)), set())

        self.assertEqual(set(mdl.get_submodels()), set(self.submodels))
        self.assertEqual(set(mdl.get_submodels(__type=Submodel)), set(self.submodels))
        self.assertEqual(set(mdl.get_submodels(__type=SpeciesType)), set())

        self.assertEqual(set(mdl.get_species()), set(self.species))
        self.assertEqual(set(mdl.get_species(__type=Species)), set(self.species))
        self.assertEqual(set(mdl.get_species(__type=Submodel)), set())

        self.assertEqual(set(mdl.get_distribution_init_concentrations()), set(self.distribution_init_concentrations))
        self.assertEqual(set(mdl.get_distribution_init_concentrations(__type=DistributionInitConcentration)),
                         set(self.distribution_init_concentrations))
        self.assertEqual(set(mdl.get_distribution_init_concentrations(__type=Submodel)), set())

        self.assertEqual(set(mdl.get_reactions()), set(self.reactions))
        self.assertEqual(set(mdl.get_reactions(__type=Reaction)), set(self.reactions))
        self.assertEqual(set(mdl.get_reactions(__type=Submodel)), set())

        self.assertEqual(set(mdl.get_rate_laws()), set(self.rate_laws))
        self.assertEqual(set(mdl.get_rate_laws(__type=RateLaw)), set(self.rate_laws))
        self.assertEqual(set(mdl.get_rate_laws(__type=Submodel)), set())

        self.assertEqual(set(mdl.get_parameters()), set(self.parameters))
        self.assertEqual(set(mdl.get_parameters(__type=Parameter)), set(self.parameters))
        self.assertEqual(set(mdl.get_parameters(__type=Submodel)), set())

        self.assertEqual(set(mdl.get_references()), set(self.references))
        self.assertEqual(set(mdl.get_references(__type=Reference)), set(self.references))
        self.assertEqual(set(mdl.get_references(__type=Submodel)), set())

        self.assertNotEqual(set(mdl.get_dfba_obj_reactions(__type=DfbaObjReaction)), set())
        self.assertEqual(set(mdl.get_dfba_obj_reactions(__type=Reaction)), set())

        self.assertEqual(set(self.dfba_obj.get_products()), set([
            self.species[3],
            self.species[4],
            self.species[1],
        ]))
        self.assertEqual(set(self.dfba_obj.get_products(__type=Species)), set([
            self.species[3],
            self.species[4],
            self.species[1],
        ]))
        self.assertEqual(set(self.dfba_obj.get_products(__type=Reaction)), set())

    def test_get_components_2(self):
        model = self.model

        self.assertEqual(model.get_components(id='comp_0'), [self.comp_0])
        self.assertEqual(model.get_components(__type=Compartment, id='comp_0'), [self.comp_0])
        self.assertEqual(model.get_components(__type=SpeciesType, id='spec_type_1'), [self.species_types[1]])
        self.assertEqual(model.get_components(__type=Submodel, id='submodel_1'), [self.submdl_1])
        self.assertEqual(model.get_components(__type=Reaction, id='rxn_1'), [self.rxn_1])
        self.assertEqual(model.get_components(__type=Parameter, id='param_2'), [self.parameters[2]])
        self.assertEqual(model.get_components(__type=Reference, id='ref_1'), [self.references[1]])
        self.assertEqual(model.get_components(__type=Reaction, id='rxn_3'), [])

    def test_validate_acyclic_compartments(self):
        model = Model(id='model', version='0.0.1')
        e = model.compartments.create(id='e')
        c = model.compartments.create(id='c', parent_compartment=e)
        m = model.compartments.create(id='m', parent_compartment=c)

        self.assertEqual(model.validate(), None)

        e.parent_compartment = m
        self.assertNotEqual(model.validate(), None)

    def test_get_root_cellular_compartments(self):
        model = Model()
        extracellular = model.compartments.create(id='extracellular', biological_type=onto['WC:extracellular_compartment'])
        cytosol = model.compartments.create(
            id='cytosol', biological_type=onto['WC:cellular_compartment'], parent_compartment=extracellular)
        nucleus = model.compartments.create(
            id='nucleus', biological_type=onto['WC:cellular_compartment'], parent_compartment=cytosol)
        nucleus_dna = model.compartments.create(
            id='nucleus_dna', biological_type=onto['WC:cellular_compartment'], parent_compartment=nucleus)
        mitochondria = model.compartments.create(
            id='mitochondria', biological_type=onto['WC:cellular_compartment'], parent_compartment=cytosol)
        mitochondria_dna = model.compartments.create(
            id='mitochondria_dna', biological_type=onto['WC:cellular_compartment'], parent_compartment=mitochondria)

        self.assertEqual(set(model.get_root_compartments()), set([extracellular, cytosol]))
        self.assertEqual(model.get_root_compartments(biological_type=onto['WC:cellular_compartment']), [cytosol])
        self.assertEqual(model.get_root_compartments(biological_type=onto['WC:extracellular_compartment']), [extracellular])

    def test_get_sub_compartments(self):
        model = Model()
        extracellular = model.compartments.create(id='extracellular', biological_type=onto['WC:extracellular_compartment'])
        cytosol = model.compartments.create(
            id='cytosol', biological_type=onto['WC:cellular_compartment'], parent_compartment=extracellular)
        nucleus = model.compartments.create(
            id='nucleus', biological_type=onto['WC:cellular_compartment'], parent_compartment=cytosol)
        nucleus_dna = model.compartments.create(
            id='nucleus_dna', biological_type=onto['WC:cellular_compartment'], parent_compartment=nucleus)
        mitochondria = model.compartments.create(
            id='mitochondria', biological_type=onto['WC:cellular_compartment'], parent_compartment=cytosol)
        mitochondria_dna = model.compartments.create(
            id='mitochondria_dna', biological_type=onto['WC:cellular_compartment'], parent_compartment=mitochondria)

        self.assertEqual(extracellular.get_sub_compartments(), [
            cytosol,
        ])
        self.assertEqual(len(extracellular.get_sub_compartments(nested=True)), 5)
        self.assertEqual(set(extracellular.get_sub_compartments(nested=True)), set([
            cytosol,
            nucleus,
            nucleus_dna,
            mitochondria,
            mitochondria_dna,
        ]))

        self.assertEqual(set(cytosol.get_sub_compartments()), set([
            nucleus,
            mitochondria,
        ]))
        self.assertEqual(len(cytosol.get_sub_compartments(nested=True)), 4)
        self.assertEqual(set(cytosol.get_sub_compartments(nested=True)), set([
            nucleus,
            nucleus_dna,
            mitochondria,
            mitochondria_dna,
        ]))

        self.assertEqual(nucleus.get_sub_compartments(), [
            nucleus_dna,
        ])
        self.assertEqual(len(nucleus.get_sub_compartments(nested=True)), 1)
        self.assertEqual(nucleus.get_sub_compartments(nested=True), [
            nucleus_dna,
        ])

    def test_get_tot_mean_init_volume(self):
        model = Model()
        extracellular = model.compartments.create(id='extracellular', biological_type=onto['WC:extracellular_compartment'],
                                                  init_volume=InitVolume(mean=1.))
        cytosol = model.compartments.create(
            id='cytosol', biological_type=onto['WC:cellular_compartment'], parent_compartment=extracellular,
            init_volume=InitVolume(mean=2.))
        nucleus = model.compartments.create(id='nucleus', biological_type=onto['WC:cellular_compartment'],
                                            parent_compartment=cytosol,
                                            init_volume=InitVolume(mean=3.))
        nucleus_dna = model.compartments.create(
            id='nucleus_dna', biological_type=onto['WC:cellular_compartment'], parent_compartment=nucleus,
            init_volume=InitVolume(mean=4.))
        mitochondria = model.compartments.create(
            id='mitochondria', biological_type=onto['WC:cellular_compartment'], parent_compartment=cytosol,
            init_volume=InitVolume(mean=5.))
        mitochondria_dna = model.compartments.create(
            id='mitochondria_dna', biological_type=onto['WC:cellular_compartment'], parent_compartment=mitochondria,
            init_volume=InitVolume(mean=6.))
        self.assertEqual(extracellular.get_tot_mean_init_volume(), 21.)
        self.assertEqual(cytosol.get_tot_mean_init_volume(), 20.)
        self.assertEqual(nucleus.get_tot_mean_init_volume(), 7.)
        self.assertEqual(nucleus_dna.get_tot_mean_init_volume(), 4.)
        self.assertEqual(mitochondria.get_tot_mean_init_volume(), 11.)
        self.assertEqual(mitochondria_dna.get_tot_mean_init_volume(), 6.)

    def test_species_type_has_carbon(self):
        'C' in self.species_types[0].structure.empirical_formula
        self.assertFalse(self.species_types[0].structure.has_carbon())
        self.assertTrue(self.species_types[1].structure.has_carbon())

    def test_species_gen_id(self):
        self.assertEqual(self.species[3].gen_id(), 'spec_type_3[comp_1]')

    def test_species__gen_id(self):
        self.assertEqual(Species._gen_id('spec_type_3', 'comp_1'), 'spec_type_3[comp_1]')

    def test_species_parse_id(self):
        self.assertEqual(Species.parse_id('spec_type_3[comp_1]'), ('spec_type_3', 'comp_1'))
        self.assertEqual(Species.parse_id('1st[comp_1]'), ('1st', 'comp_1'))
        self.assertEqual(Species.parse_id('1st[1comp]'), ('1st', '1comp'))
        self.assertEqual(Species.parse_id('1ST[1Comp]'), ('1ST', '1Comp'))
        with self.assertRaisesRegex(ValueError, ''):
            Species.parse_id('123[comp_1]')
        with self.assertRaisesRegex(ValueError, ''):
            Species.parse_id('spec_type_3[123]')
        with self.assertRaisesRegex(ValueError, ''):
            Species.parse_id('spec_type_3')

    def test_species_get(self):
        self.assertEqual(Species.get([], self.species), [])
        self.assertEqual(Species.get(['X'], self.species), [None])
        self.assertEqual(Species.get(['spec_type_0[comp_0]'], self.species), [self.species[0]])
        ids = ["spec_type_{}[comp_0]".format(i) for i in range(4, 8)]
        self.assertEqual(Species.get(ids, self.species), self.species[4:])
        ids.append('X')
        self.assertEqual(Species.get(ids, self.species), self.species[4:] + [None])

    def test_distribution_init_concentration_serialize(self):
        self.assertEqual(self.distribution_init_concentrations[0].serialize(), 'dist-init-conc-spec_type_0[comp_0]')
        self.assertEqual(self.distribution_init_concentrations[1].serialize(), 'dist-init-conc-spec_type_1[comp_0]')
        self.assertEqual(self.distribution_init_concentrations[2].serialize(), 'dist-init-conc-spec_type_2[comp_0]')
        self.assertEqual(self.distribution_init_concentrations[3].serialize(), 'dist-init-conc-spec_type_3[comp_1]')

    def test_distribution_init_concentration_validate(self):
        conc = DistributionInitConcentration(id='dist-init-conc-species_0', species=Species(id='species_0'), mean=1.)
        self.assertEqual(conc.validate(), None)

        conc = DistributionInitConcentration(id='dist-init-conc-species_0', species=Species(id='species_0'), mean=-1.)
        self.assertNotEqual(conc.validate(), None)

        conc = DistributionInitConcentration(id='dist-init-conc-species-0', species=Species(id='species_0'), mean=1.)
        self.assertNotEqual(conc.validate(), None)

    def test_reaction_participant_serialize(self):
        self.assertEqual(set([part.serialize() for part in self.rxn_0.participants]), set([
            '(-2) spec_type_0[comp_0]', '(-3.500000e+00) spec_type_1[comp_0]', 'spec_type_2[comp_0]'
        ]))

    def test_validate_reaction_balance(self):
        c = Compartment()
        st_1 = SpeciesType(structure=ChemicalStructure(empirical_formula=EmpiricalFormula('CH1N2OP2'), charge=1))
        st_2 = SpeciesType(structure=ChemicalStructure(empirical_formula=EmpiricalFormula('C2H2N4O2P4'), charge=2))
        st_3 = SpeciesType(structure=ChemicalStructure(empirical_formula=EmpiricalFormula('C3H3N6O3P6'), charge=3))
        st_4 = SpeciesType(structure=ChemicalStructure(empirical_formula=EmpiricalFormula('CH1N2'), charge=2))
        st_5 = SpeciesType(structure=ChemicalStructure(empirical_formula=EmpiricalFormula('OP2'), charge=-1))
        s_1 = Species(species_type=st_1, compartment=c)
        s_2 = Species(species_type=st_2, compartment=c)
        s_3 = Species(species_type=st_3, compartment=c)
        s_4 = Species(species_type=st_4, compartment=c)
        s_5 = Species(species_type=st_5, compartment=c)

        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_lang__DOT__validation__DOT__validate_element_charge_balance', '1')
        with env:
            rxn = Reaction(id='rxn')
            rxn.participants.create(species=s_1, coefficient=-2.)
            rxn.participants.create(species=s_2, coefficient=1.)
            rv = rxn.validate()
            self.assertEqual(rv, None, str(rv))

            rxn = Reaction(id='rxn')
            rxn.participants.create(species=s_1, coefficient=-3.)
            rxn.participants.create(species=s_3, coefficient=1.)
            rv = rxn.validate()
            self.assertEqual(rv, None, str(rv))

            rxn = Reaction(id='rxn')
            rxn.participants.create(species=s_1, coefficient=-1.)
            rxn.participants.create(species=s_2, coefficient=-1.)
            rxn.participants.create(species=s_3, coefficient=1.)
            rv = rxn.validate()
            self.assertEqual(rv, None, str(rv))

            rxn = Reaction(id='rxn')
            rxn.participants.create(species=s_1, coefficient=-1.)
            rxn.participants.create(species=s_4, coefficient=1.)
            rxn.participants.create(species=s_5, coefficient=1.)
            rv = rxn.validate()
            self.assertEqual(rv, None, str(rv))

            rxn = Reaction(id='rxn')
            rxn.participants.create(species=s_4, coefficient=1.)
            rxn.participants.create(species=s_5, coefficient=1.)
            rv = rxn.validate()
            self.assertRegex(str(rv), 'element imbalanced')
            self.assertRegex(str(rv), 'charge imbalanced')

            rxn = Reaction(id='rxn')
            rxn.participants.create(species=s_1, coefficient=-3.3333333333)
            rxn.participants.create(species=s_3, coefficient=1.11111111111)
            rv = rxn.validate()
            self.assertEqual(rv, None, str(rv))

            st_1.structure.empirical_formula = EmpiricalFormula('CH1N2OP2')
            st_1.structure.charge = None
            rxn = Reaction(id='rxn')
            rxn.participants.create(species=s_1, coefficient=-2.)
            rxn.participants.create(species=s_2, coefficient=1.)
            rv = rxn.validate()
            self.assertRegex(str(rv), 'Charge must be defined ')

            with self.assertRaisesRegex(ValueError, 'not a valid formula'):
                st_1.structure.empirical_formula = EmpiricalFormula('1')

        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_lang__DOT__validation__DOT__validate_element_charge_balance', '0')
        with env:
            rv = rxn.validate()
            self.assertEqual(rv, None, str(rv))

    def test_reaction_validate(self):
        c = Compartment()
        d = Compartment()
        st = SpeciesType(structure=ChemicalStructure(empirical_formula=EmpiricalFormula('CHO'), charge=1))
        spec_c = Species(species_type=st, compartment=c)
        spec_d = Species(species_type=st, compartment=d)
        rxn = Reaction(id='rxn', reversible=True,
                       flux_bounds=FluxBounds(min=-1., max=1., units=unit_registry.parse_units('M s^-1')),
                       participants=[
                           SpeciesCoefficient(species=spec_c, coefficient=-1.),
                           SpeciesCoefficient(species=spec_d, coefficient=1.),
                       ],
                       rate_laws=[
                           RateLaw(direction=RateLawDirection.forward),
                           RateLaw(direction=RateLawDirection.backward),
                       ])
        rv = rxn.validate()
        self.assertEqual(rv, None, str(rv))

        rxn.flux_bounds.units = None
        rv = rxn.validate()
        self.assertNotEqual(rv, None, str(rv))

    def test_rate_gen_id(self):
        self.assertEqual(self.rate_laws[0].id, 'rxn_0-forward')
        self.assertEqual(self.rate_laws[1].id, 'rxn_1-forward')
        self.assertEqual(self.rate_laws[2].id, 'rxn_2-forward')

    def test_rate_law_expression_serialize(self):
        self.assertEqual(self.rate_laws[0].expression.serialize(),
                         'k_cat_0 * {0} / (k_m_0 + {0})'.format(self.species[5].serialize()))
        self.assertEqual(self.rate_laws[1].expression.serialize(),
                         'k_cat_1 * {0} / (k_m_1 + {0})'.format(self.species[6].serialize()))
        self.assertEqual(self.rate_laws[2].expression.serialize(),
                         '{1} * {0} / ({2} + {0})'.format(
            self.species[7].serialize(),
            self.parameters[0].id,
            self.parameters[1].id,
        ))

    def test_rate_law_expression_deserialize(self):
        objs = {
            SpeciesType: {
                'spec_0': SpeciesType(id='spec_0'),
                'spec_1': SpeciesType(id='spec_1'),
                'spec_2': SpeciesType(id='spec_2')},
            Compartment: {
                'c_0': Compartment(id='c_0'),
                'c_1': Compartment(id='c_1'),
                'c_2': Compartment(id='c_2')
            },
            Parameter: {
                'p_1': Parameter(id='p_1'),
                'p_2': Parameter(id='p_2'),
            },
            Species: {
            },
            Parameter: {
                'k_cat': Parameter(id='k_cat', value=1),
                'k_m': Parameter(id='k_m', value=2),
            }
        }
        objs[Species]['spec_0[c_0]'] = Species(id='spec_0[c_0]',
                                               species_type=objs[SpeciesType]['spec_0'],
                                               compartment=objs[Compartment]['c_0'])
        objs[Species]['spec_0[c_1]'] = Species(id='spec_0[c_1]',
                                               species_type=objs[SpeciesType]['spec_0'],
                                               compartment=objs[Compartment]['c_1'])
        objs[Species]['spec_2[c_1]'] = Species(id='spec_2[c_1]',
                                               species_type=objs[SpeciesType]['spec_2'],
                                               compartment=objs[Compartment]['c_1'])
        objs[Species]['spec_1[c_1]'] = Species(id='spec_1[c_1]',
                                               species_type=objs[SpeciesType]['spec_1'],
                                               compartment=objs[Compartment]['c_1'])

        expression = 'spec_0[c_0]'
        expression1, error = RateLawExpression.deserialize(expression, objs)
        self.assertEqual(error, None)
        self.assertEqual(expression1.expression, expression)
        self.assertEqual(expression1.species, [objs[Species]['spec_0[c_0]']])
        self.assertEqual(expression1.parameters, [])
        self.assertEqual(set(objs[RateLawExpression].values()), set([expression1]))

        expression = 'spec_0[c_1] / (spec_2[c_1])'
        expression2, error = RateLawExpression.deserialize(expression, objs)
        self.assertEqual(error, None)
        self.assertEqual(expression2.expression, expression)
        self.assertEqual(set(expression2.species), set([objs[Species]['spec_0[c_1]'], objs[Species]['spec_2[c_1]']]))
        self.assertEqual(expression2.parameters, [])
        self.assertEqual(set(objs[RateLawExpression].values()), set([expression1, expression2]))

        expression = 'spec_0[c_3] / (spec_1[c_1])'
        expression, error = RateLawExpression.deserialize(expression, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(expression, None)
        self.assertEqual(set(objs[RateLawExpression].values()), set([expression1, expression2]))

        expression = 'spec_3[c_0] / (spec_1[c_1])'
        expression, error = RateLawExpression.deserialize(expression, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(expression, None)
        self.assertEqual(set(objs[RateLawExpression].values()), set([expression1, expression2]))

        # exception
        expression3, error = RateLawExpression.deserialize('2', objs)
        self.assertEqual(error, None)
        self.assertEqual(expression3.expression, '2')

        # with parameters
        expression = 'k_cat * spec_0[c_1] / (k_m + spec_2[c_1])'
        expression4, error = RateLawExpression.deserialize(expression, objs)
        self.assertEqual(error, None)
        self.assertEqual(expression4.expression, expression)
        self.assertEqual(set(expression4.species), set([objs[Species]['spec_0[c_1]'], objs[Species]['spec_2[c_1]']]))
        self.assertEqual(set(expression4.parameters), set(objs[Parameter].values()))
        self.assertEqual(set(objs[RateLawExpression].values()), set([expression1, expression2, expression3, expression4]))

        expression = 'p_1 * spec_0[c_1] / (p_2 + spec_2[c_1])'
        expression, error = RateLawExpression.deserialize(expression, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(expression, None)
        self.assertEqual(set(objs[RateLawExpression].values()), set([expression1, expression2, expression3, expression4]))

    def test_rate_law_validate(self):
        species_types = [
            SpeciesType(id='spec_0'),
            SpeciesType(id='spec_1'),
        ]
        compartments = [
            Compartment(id='c_0'),
            Compartment(id='c_1'),
        ]
        parameters = [
            Parameter(id='p_0', value=1.),
        ]

        # unknown specie error
        expression = 'spec_x[c_0]'
        expression = RateLawExpression(
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0])
            ])
        rate_law = RateLaw(
            id='rxn-forward',
            reaction=Reaction(id='rxn'),
            expression=expression,
        )
        self.assertNotEqual(rate_law.expression.validate(), None)

        # unknown parameter error
        expression = 'p_x'
        expression = RateLawExpression(
            expression=expression,
            parameters=[parameters[0]])
        rate_law = RateLaw(
            id='rxn-forward',
            reaction=Reaction(id='rxn'),
            expression=expression,
        )
        self.assertNotEqual(rate_law.expression.validate(), None)

        # Name error
        expression = 'not_k_cat * spec_0[c_0]'
        expression = RateLawExpression(
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0])
            ])
        rate_law = RateLaw(
            id='rxn-forward',
            reaction=Reaction(id='rxn'),
            expression=expression
        )
        self.assertNotEqual(rate_law.expression.validate(), None)

        # syntax error
        expression = '* spec_0[c_0]'
        expression = RateLawExpression(
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0])
            ])
        rate_law = RateLaw(
            id='rxn-forward',
            reaction=Reaction(id='rxn'),
            expression=expression
        )
        self.assertRegex(str(rate_law.expression.validate()), 'SyntaxError')

        # No error
        expression = 'k_cat * spec_0[c_0]'
        expression, errors = RateLawExpression.deserialize(expression, {
            Species: {
                'spec_0[c_0]': Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0]),
            },
            Parameter: {
                'k_cat': Parameter(id='k_cat', value=1, units=unit_registry.parse_units('molecule^-1 s^-1')),
            },
        })
        rate_law = RateLaw(
            id='rxn-forward',
            reaction=Reaction(id='rxn'),
            expression=expression,
            units=unit_registry.parse_units('s^-1'),
        )
        error = rate_law.validate()
        self.assertEqual(error, None, str(error))
        error = rate_law.expression.validate()
        self.assertEqual(error, None, str(error))

        # No error with parameters
        expression = 'p_0 * spec_0[c_0]'
        expression, _ = RateLawExpression.deserialize(expression, {
            Species: {'spec_0[c_0]': Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0])},
            Parameter: {'p_0': Parameter(id='p_0', value=1, units=unit_registry.parse_units('molecule^-1 s^-1'))},
        })
        rate_law = RateLaw(
            id='rxn-forward',
            reaction=Reaction(id='rxn'),
            expression=expression,
            units=unit_registry.parse_units('s^-1'),
        )
        error = rate_law.validate()
        self.assertEqual(error, None, str(error))
        error = rate_law.expression.validate()
        self.assertEqual(error, None, str(error))

        # positive example
        expression, _ = RateLawExpression.deserialize('p', {
            Parameter: {'p': Parameter(id='p', value=1., units=unit_registry.parse_units('s^-1'))},
        })
        rate_law = RateLaw(
            id='rxn-forward',
            reaction=Reaction(id='rxn'),
            expression=expression,
            units=unit_registry.parse_units('s^-1'),
        )
        error = rate_law.validate()
        self.assertEqual(error, None, str(error))

        # invalid id
        rate_law.id = 'a' * 1000
        error = rate_law.validate()
        self.assertNotEqual(error, None, str(error))

        rate_law.id = 'rxn-backward'
        error = rate_law.validate()
        self.assertNotEqual(error, None, str(error))

        rate_law.id = 'rxn-forward'
        rate_law.expression._parsed_expression._compiled_expression_with_units = 'p'
        error = rate_law.validate()
        self.assertNotEqual(error, None, str(error))

        rate_law.id = 'rxn-forward'
        rate_law.expression._parsed_expression._compiled_expression_with_units = 'True'
        error = rate_law.validate()
        self.assertNotEqual(error, None, str(error))

        rate_law.id = 'rxn-forward'
        rate_law.expression = None
        error = rate_law.validate()
        self.assertNotEqual(error, None, str(error))

        # no parsed expression
        expression, _ = RateLawExpression.deserialize('p', {
            Parameter: {'p': Parameter(id='p', value=1., units=unit_registry.parse_units('s^-1'))},
        })
        expression._parsed_expression = None
        rate_law = RateLaw(
            id='rxn-forward',
            reaction=Reaction(id='rxn'),
            expression=expression,
            units=unit_registry.parse_units('s^-1'),
        )
        error = rate_law.validate()
        self.assertNotEqual(error, None, str(error))

    def test_rate_law_expression_validate(self):
        species_types = [
            SpeciesType(id='spec_0'),
            SpeciesType(id='spec_1'),
            SpeciesType(id='spec_2'),
        ]
        compartments = [
            Compartment(id='c_0'),
            Compartment(id='c_1'),
            Compartment(id='c_2'),
        ]
        parameters = [
            Parameter(id='p_0', value=1., units=unit_registry.parse_units('molecule^-1 s^-1')),
            Parameter(id='p_1', value=1.),
            Parameter(id='k_m', value=1.),
        ]

        expression = 'spec_0[c_0]'
        expression = RateLawExpression(
            rate_laws=[RateLaw()],
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0])
            ])
        error = expression.validate()
        self.assertEqual(error, None, str(error))

        expression = 'spec_0[c_0] * spec_1[c_2]'
        expression = RateLawExpression(
            rate_laws=[RateLaw()],
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0]),
                Species(id='spec_1[c_2]', species_type=species_types[1], compartment=compartments[2]),
            ])
        error = expression.validate()
        self.assertEqual(error, None, str(error))

        expression = 'spec_0[c_0] * spec_1[c_2]'
        expression = RateLawExpression(
            rate_laws=[RateLaw()],
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0]),
                Species(id='spec_1[c_1]', species_type=species_types[1], compartment=compartments[1]),
                Species(id='spec_1[c_2]', species_type=species_types[1], compartment=compartments[2]),
            ])
        self.assertNotEqual(expression.validate(), None)

        expression = 'spec_0[c_0] * spec_1[c_2]'
        expression = RateLawExpression(
            rate_laws=[RateLaw()],
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0]),
            ])
        self.assertNotEqual(expression.validate(), None)

        # parameters
        expression = 'p_0 * spec_0[c_0]'
        expression = RateLawExpression(
            rate_laws=[RateLaw()],
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0])
            ],
            parameters=[parameters[0]])
        error = expression.validate()
        self.assertEqual(error, None, str(error))

        expression = 'p_0 * spec_0[c_0]'
        expression = RateLawExpression(
            rate_laws=[RateLaw()],
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0])
            ],
            parameters=[parameters[0], parameters[1]])
        self.assertNotEqual(expression.validate(), None)

        expression = 'p_1 * spec_0[c_0]'
        expression = RateLawExpression(
            rate_laws=[RateLaw()],
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0])
            ],
            parameters=[])
        self.assertNotEqual(expression.validate(), None)

        expression = 'k_m * spec_0[c_0]'
        expression = RateLawExpression(
            rate_laws=[RateLaw()],
            expression=expression,
            species=[
                Species(id='spec_0[c_0]', species_type=species_types[0], compartment=compartments[0])
            ],
            parameters=[parameters[2]])
        error = expression.validate()
        self.assertEqual(error, None, str(error))

    def test_rate_law_species(self):
        self.assertEqual(self.rxn_0.rate_laws[0].expression.species, self.species[5:6])
        self.assertEqual(self.rxn_1.rate_laws[0].expression.species, self.species[6:7])
        self.assertEqual(self.rxn_2.rate_laws[0].expression.species, self.species[7:8])

    def test_rate_law_parameters(self):
        self.assertEqual([p.id for p in self.rxn_0.rate_laws[0].expression.parameters], ['k_cat_0', 'k_m_0'])
        self.assertEqual([p.id for p in self.rxn_1.rate_laws[0].expression.parameters], ['k_cat_1', 'k_m_1'])
        self.assertEqual(self.rxn_2.rate_laws[0].expression.parameters, self.parameters[0:2])

    def test_parameter_validate_unique(self):
        self.assertEqual(Parameter.validate_unique(self.parameters), None)

        model = Model()
        params = [
            Parameter(id='a', model=model),
            Parameter(id='b', model=model),
        ]
        self.assertEqual(Parameter.validate_unique(params), None)

        model = Model()
        params = [
            Parameter(id='a', model=model),
            Parameter(id='a', model=model),
        ]
        self.assertNotEqual(Parameter.validate_unique(params), None)

        submodel = Submodel()
        params = [
            Parameter(id='a'),
            Parameter(id='b'),
        ]
        self.assertEqual(Parameter.validate_unique(params), None)

        submodel = Submodel()
        params = [
            Parameter(id='a'),
            Parameter(id='a'),
        ]
        self.assertNotEqual(Parameter.validate_unique(params), None)

        model = Model()
        submodel = Submodel()
        params = [
            Parameter(id='a', model=model),
            Parameter(id='a', model=model),
        ]
        self.assertNotEqual(Parameter.validate_unique(params), None)

        params = [
            Parameter(id='a'),
            Parameter(id='b'),
        ]
        self.assertEqual(Parameter.validate_unique(params), None)

        params = [
            Parameter(id='a', model=model),
            Parameter(id='b', model=model),
        ]
        self.assertEqual(Parameter.validate_unique(params), None)

    def test_identifier_serialize(self):
        self.assertEqual(self.identifiers[0].serialize(), '{}: {}'.format('x', 'y'))
        self.assertEqual(self.identifiers[1].serialize(), '{}: {}'.format('x', 'yy'))
        self.assertEqual(self.identifiers[2].serialize(), '{}: {}'.format('x', 'yyy'))

    def test_identifier_deserialize(self):
        objs = {
            Identifier: {
            }
        }
        identifiers, errors = Model.Meta.attributes['identifiers'].deserialize('db_1: id_1, db_2: id_2', objs)
        self.assertEqual(len(identifiers), 2)
        self.assertEqual(errors, None)
        self.assertEqual(len(objs[Identifier]), 2)

        identifiers, errors = Submodel.Meta.attributes['identifiers'].deserialize('db_2: id_2, db_3: id_3, db_4: id_4', objs)
        self.assertEqual(len(identifiers), 3)
        self.assertEqual(errors, None)
        self.assertEqual(len(objs[Identifier]), 4)

        identifiers, errors = Model.Meta.attributes['identifiers'].deserialize('db_1', objs)
        self.assertNotEqual(errors, None)

        identifiers, errors = Submodel.Meta.attributes['identifiers'].deserialize('db_1', objs)
        self.assertNotEqual(errors, None)

    def test_ReactionParticipantAttribute_serialize(self):
        attr = ReactionParticipantAttribute()
        self.assertEqual(attr.serialize(self.rxn_0.participants),
                         '[comp_0]: (2) spec_type_0 + (3.500000e+00) spec_type_1 ==> spec_type_2')
        self.assertEqual(attr.serialize(self.rxn_1.participants),
                         '(2) spec_type_0[comp_0] + (3) spec_type_1[comp_0] ==> (2) spec_type_3[comp_1]')
        self.assertEqual(attr.serialize(None), '')

    def test_ReactionParticipantAttribute_deserialize(self):
        objs = {
            SpeciesType: {
                'spec_0': SpeciesType(id='spec_0'),
                'spec_1': SpeciesType(id='spec_1'),
                'spec_2': SpeciesType(id='spec_2')},
            Compartment: {
                'c_0': Compartment(id='c_0'),
                'c_1': Compartment(id='c_1'),
                'c_2': Compartment(id='c_2')
            },
            Species: {
            },
        }
        objs[Species]['spec_0[c_0]'] = Species(
            species_type=objs[SpeciesType]['spec_0'],
            compartment=objs[Compartment]['c_0'])
        objs[Species]['spec_1[c_0]'] = Species(
            species_type=objs[SpeciesType]['spec_1'],
            compartment=objs[Compartment]['c_0'])
        objs[Species]['spec_2[c_0]'] = Species(
            species_type=objs[SpeciesType]['spec_2'],
            compartment=objs[Compartment]['c_0'])
        objs[Species]['spec_2[c_1]'] = Species(
            species_type=objs[SpeciesType]['spec_2'],
            compartment=objs[Compartment]['c_1'])
        for species in objs[Species].values():
            species.id = species.gen_id()

        attr = ReactionParticipantAttribute()

        parts1, error = attr.deserialize('[c_0]: (2) spec_0 + (3.5) spec_1 ==> spec_2', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts1]), set([
            '(-2) spec_0[c_0]',
            '(-3.500000e+00) spec_1[c_0]',
            'spec_2[c_0]',
        ]))
        self.assertEqual(len(objs[SpeciesCoefficient]), 3)
        self.assertEqual(set(objs[SpeciesCoefficient].values()), set(parts1))
        self.assertEqual(len(objs[Species]), 4)

        parts2, error = attr.deserialize(
            '(2) spec_0[c_0] + (3) spec_1[c_0] ==> (2) spec_2[c_1]', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts2]), set(
            ['(-2) spec_0[c_0]', '(-3) spec_1[c_0]', '(2) spec_2[c_1]']))
        self.assertEqual(set([p.serialize() for p in objs[SpeciesCoefficient].values()]),
                         set([p.serialize() for p in parts1 + parts2]))
        self.assertEqual(len(objs[Species]), 4)

        # negative examples
        parts3, error = attr.deserialize(
            '(2) spec_0[c_0] + (3) spec_1[c_0] ==> (2) spec_2[c_3]', objs)
        self.assertNotEqual(error, None)
        self.assertEqual(parts3, None)

        parts3, error = attr.deserialize(
            '(2) spec_0[c_0] + (3) spec_1[c_0] => (2) spec_2[c_1]', objs)
        self.assertNotEqual(error, None)
        self.assertEqual(parts3, None)

        parts3, error = attr.deserialize(
            '(2) spec_0[c_0] + (3) spec_1[c_0] ==> (-2) spec_2[c_1]', objs)
        self.assertNotEqual(error, None)
        self.assertEqual(parts3, None)

        parts3, error = attr.deserialize(
            '[c_0]: (2) spec_0[c_0] + (3) spec_1[c_0] ==> (2) spec_2[c_1]', objs)
        self.assertNotEqual(error, None)
        self.assertEqual(parts3, None)

        parts3, error = attr.deserialize(
            '[c_0]: (2) spec_0 + (3) spec_1 ==> (2) spec_3', objs)
        self.assertNotEqual(error, None)
        self.assertEqual(parts3, None)

        parts, error = attr.deserialize('[c_3]: (2) spec_0 + (3.5) spec_1 ==> spec_2', objs)
        self.assertNotEqual(error, None)
        self.assertEqual(parts, None)

        # empty LHS
        parts, error = attr.deserialize('==> spec_2[c_1]', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts]), set(
            ['spec_2[c_1]']))

        parts, error = attr.deserialize('[c_1]: ==> spec_2', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts]), set(
            ['spec_2[c_1]']))

        parts, error = attr.deserialize('[c_1]:  ==>spec_2 ', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts]), set(
            ['spec_2[c_1]']))

        # empty RHS
        parts, error = attr.deserialize('spec_2[c_1] ==>', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts]), set(
            ['(-1) spec_2[c_1]']))

        parts, error = attr.deserialize('[c_1]: spec_2 ==>', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts]), set(
            ['(-1) spec_2[c_1]']))

        parts, error = attr.deserialize('[c_1]:spec_2==> ', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts]), set(
            ['(-1) spec_2[c_1]']))

        # both empty
        parts, error = attr.deserialize('==>', objs)
        self.assertEqual(error, None)
        self.assertEqual(parts, [])

        parts, error = attr.deserialize('[c_1]: ==>', objs)
        self.assertEqual(error, None)
        self.assertEqual(parts, [])

        parts, error = attr.deserialize('[c_1]:  ==>  ', objs)
        self.assertEqual(error, None)
        self.assertEqual(parts, [])

        parts, error = attr.deserialize('[c_1]:==>', objs)
        self.assertEqual(error, None)
        self.assertEqual(parts, [])

        # repeated species
        objs[Species] = {
            'spec_2[c_1]': Species(id='spec_2[c_1]', species_type=objs[SpeciesType]['spec_2'], compartment=objs[Compartment]['c_1']),
        }
        objs[SpeciesCoefficient] = {
            '(-1) spec_2[c_1]': SpeciesCoefficient(species=objs[Species]['spec_2[c_1]'], coefficient=-1),
            'spec_2[c_1]': SpeciesCoefficient(species=objs[Species]['spec_2[c_1]'], coefficient=1),
        }
        parts, error = attr.deserialize('[c_1]: spec_2 ==> spec_2 + spec_2', objs)
        self.assertNotEqual(error, None)
        self.assertEqual(parts, None)

        parts, error = attr.deserialize('spec_2[c_1] ==> spec_2[c_1] + spec_2[c_1]', objs)
        self.assertNotEqual(error, None)
        self.assertEqual(parts, None)

        # species type whose id starts with number
        objs = {
            SpeciesType: {'1st': SpeciesType(id='1st'), '2st': SpeciesType(id='2st')},
            Compartment: {'1comp': Compartment(id='1comp'), '2comp': Compartment(id='2comp')},
            Species: {},
        }
        objs[Species]['1st[1comp]'] = Species(
            id='1st[1comp]',
            species_type=objs[SpeciesType]['1st'],
            compartment=objs[Compartment]['1comp'])
        objs[Species]['1st[2comp]'] = Species(
            id='1st[2comp]',
            species_type=objs[SpeciesType]['1st'],
            compartment=objs[Compartment]['2comp'])
        objs[Species]['2st[1comp]'] = Species(
            id='2st[1comp]',
            species_type=objs[SpeciesType]['2st'],
            compartment=objs[Compartment]['1comp'])
        objs[SpeciesCoefficient] = {
            '(-1) 1st[1comp]': SpeciesCoefficient(species=objs[Species]['1st[1comp]'], coefficient=-1),
            '1st[2comp]': SpeciesCoefficient(species=objs[Species]['1st[2comp]'], coefficient=1),
            '2st[1comp]': SpeciesCoefficient(species=objs[Species]['2st[1comp]'], coefficient=1),
        }

        parts, error = attr.deserialize('[1comp]: 1st ==> 2st', objs)
        self.assertEqual(error, None)
        self.assertEqual(set(parts), set([objs[SpeciesCoefficient]['(-1) 1st[1comp]'], objs[SpeciesCoefficient]['2st[1comp]']]))

        parts, error = attr.deserialize('1st[1comp] ==> 1st[2comp]', objs)
        self.assertEqual(error, None)
        self.assertEqual(set(parts), set([objs[SpeciesCoefficient]['(-1) 1st[1comp]'], objs[SpeciesCoefficient]['1st[2comp]']]))

    def test_ReactionParticipantAttribute_validate(self):
        species_types = [
            SpeciesType(id='A', structure=ChemicalStructure(empirical_formula=EmpiricalFormula('CHO'), charge=2)),
            SpeciesType(id='B', structure=ChemicalStructure(empirical_formula=EmpiricalFormula('C1H1O1'), charge=2)),
        ]
        compartments = [
            Compartment(id='c'),
            Compartment(id='e'),
        ]
        species = [
            Species(id='A[c]', species_type=species_types[0], compartment=compartments[0]),
            Species(id='A[e]', species_type=species_types[0], compartment=compartments[1]),
            Species(id='B[c]', species_type=species_types[1], compartment=compartments[0]),
            Species(id='B[e]', species_type=species_types[1], compartment=compartments[1]),
        ]

        attr = ReactionParticipantAttribute()
        attr.related_class = SpeciesCoefficient

        rxn = Reaction(participants=[
            SpeciesCoefficient(species=species[0], coefficient=-1),
            SpeciesCoefficient(species=species[1], coefficient=1),
        ])
        self.assertEqual(attr.validate(None, rxn.participants), None)

        self.assertNotEqual(attr.validate(None, [
            SpeciesCoefficient(species=species[0], coefficient=-1),
            SpeciesCoefficient(species=species[0], coefficient=1),
        ]), None)

        self.assertNotEqual(attr.validate(None, [
            SpeciesCoefficient(species=species[0], coefficient=-2),
            SpeciesCoefficient(species=species[0], coefficient=1),
            SpeciesCoefficient(species=species[0], coefficient=1),
        ]), None)

        self.assertNotEqual(attr.validate(None, [
            SpeciesCoefficient(species=species[0], coefficient=-2),
            SpeciesCoefficient(species=species[1], coefficient=-1),
            SpeciesCoefficient(species=species[0], coefficient=2),
            SpeciesCoefficient(species=species[1], coefficient=1),
        ]), None)

        self.assertNotEqual(attr.validate(None, []), None)

    def test_ExpressionManyToOneAttribute_serialize(self):
        rxn = self.rxn_0
        rate_law = rxn.rate_laws[0]
        expression = rate_law.expression

        attr = ExpressionManyToOneAttribute(RateLawExpression)
        self.assertEqual(attr.serialize(expression), expression.expression)

    def test_RateLawExpressionAttribute_deserialize(self):
        objs = {
            SpeciesType: {
                'spec_0': SpeciesType(id='spec_0'),
                'spec_1': SpeciesType(id='spec_1'),
                'spec_2': SpeciesType(id='spec_2')},
            Compartment: {
                'c_0': Compartment(id='c_0'),
                'c_1': Compartment(id='c_1'),
                'c_2': Compartment(id='c_2')
            },
            Species: {
            },
            Parameter: {
                'k_cat': Parameter(id='k_cat'),
            }
        }
        objs[Species]['spec_0[c_0]'] = Species(
            id='spec_0[c_0]',
            species_type=objs[SpeciesType]['spec_0'],
            compartment=objs[Compartment]['c_0'])

        expression = 'k_cat * spec_0[c_0]'
        attr = ExpressionManyToOneAttribute(RateLawExpression)
        expression1, error = attr.deserialize(expression, objs)
        self.assertEqual(error, None)
        self.assertEqual(expression1.expression, expression)
        self.assertEqual(expression1.species, [objs[Species]['spec_0[c_0]']])
        self.assertEqual(expression1.parameters, [objs[Parameter]['k_cat']])
        self.assertEqual(list(objs[RateLawExpression].values()), [expression1])

    def test_dfba_obj_deserialize(self):
        objs = {
            Reaction: {
                'reaction_0': Reaction(id='reaction_0'),
                'reaction_1': Reaction(id='reaction_1'),
                'reaction_2': Reaction(id='reaction_2'),
            },
            DfbaObjReaction: {
                'dfba_obj_reaction_0': DfbaObjReaction(id='dfba_obj_reaction_0'),
                'dfba_obj_reaction_1': DfbaObjReaction(id='dfba_obj_reaction_1'),
            },
        }

        value = None
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertEqual(of_expr, None)
        self.assertNotEqual(invalid_attribute, None)

        value = ''
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertEqual(of_expr, None)
        self.assertNotEqual(invalid_attribute, None)

        value = "2*dfba_obj_reaction_1 - reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertEqual(invalid_attribute, None)
        self.assertEqual(of_expr.reactions, [objs[Reaction]['reaction_1']])
        self.assertEqual(of_expr.dfba_obj_reactions, [objs[DfbaObjReaction]['dfba_obj_reaction_1']])
        self.assertEqual(of_expr._parsed_expression.is_linear, True)

        value = "2*dfba_obj_reaction_1 - pow( reaction_1, 1)"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertNotEqual(invalid_attribute, None)

        value = "2*dfba_obj_reaction_1 - pow( reaction_1, 2)"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertNotEqual(invalid_attribute, None)

        objs[Reaction]['dfba_obj_reaction_1'] = Reaction(id='dfba_obj_reaction_1')
        value = "2*dfba_obj_reaction_1 - pow( reaction_1, 2)"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertTrue(of_expr is None)
        self.assertIn(("contains multiple model object id matches: "
                       "'dfba_obj_reaction_1' as a DfbaObjReaction id, "
                       "'dfba_obj_reaction_1' as a Reaction id"),
                      invalid_attribute.messages[0])

        del objs[Reaction]['dfba_obj_reaction_1']
        value = "2*dfba_obj_reaction_1 - reaction_x"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertTrue(of_expr is None)
        self.assertIn("contains the identifier(s) 'reaction_x', which aren't the id(s) of an object",
                      invalid_attribute.messages[0])

        value = "2*dfba_obj_reaction_1 - pow( dfba_obj_reaction_1, 2)"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertNotEqual(invalid_attribute, None)

        value = 'dfba_obj_reaction_1 + reaction_1 + 2.0 * reaction_2'
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertEqual(invalid_attribute, None)
        self.assertTrue(of_expr._parsed_expression.is_linear)
        self.assertEqual(of_expr._parsed_expression.lin_coeffs, {
            DfbaObjReaction: {
                objs[DfbaObjReaction]['dfba_obj_reaction_1']: 1.0,
            },
            Reaction: {
                objs[Reaction]['reaction_1']: 1.0,
                objs[Reaction]['reaction_2']: 2.0,
            },
        })

        value = 'dfba_obj_reaction_1 + reaction_1 + 2.0 * reaction_2 + reaction_2 - 0.5 * reaction_2'
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertEqual(invalid_attribute, None)
        self.assertTrue(of_expr._parsed_expression.is_linear)
        self.assertEqual(of_expr._parsed_expression.lin_coeffs, {
            DfbaObjReaction: {
                objs[DfbaObjReaction]['dfba_obj_reaction_1']: 1.0,
            },
            Reaction: {
                objs[Reaction]['reaction_1']: 1.0,
                objs[Reaction]['reaction_2']: 2.5,
            },
        })

        value = 'dfba_obj_reaction_1 + reaction_1 + reaction_2 * 2.0'
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertEqual(invalid_attribute, None)
        self.assertFalse(of_expr._parsed_expression.is_linear)

        value = 'dfba_obj_reaction_1 * reaction_1'
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertEqual(invalid_attribute, None)
        self.assertFalse(of_expr._parsed_expression.is_linear)

        value = 'reaction_1 * reaction_1'
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertEqual(invalid_attribute, None)
        self.assertFalse(of_expr._parsed_expression.is_linear)

    def test_dfba_obj_deserialize_invalid_ids(self):

        objs = {
            Reaction: {
                'rxn': Reaction(id='rxn'),
            },
            DfbaObjReaction: {
                'rxn': DfbaObjReaction(id='rxn'),
            },
        }

        value = "2*rxn"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertEqual(of_expr, None)
        self.assertIn("multiple model object id matches",
                      invalid_attribute.messages[0])

    def test_dfba_obj_validate(self):
        objs = {
            DfbaObjReaction: {
                'dfba_obj_reaction_0': DfbaObjReaction(id='dfba_obj_reaction_0'),
                'dfba_obj_reaction_1': DfbaObjReaction(id='dfba_obj_reaction_1'),
            },
            Reaction: {
                'reaction_0': Reaction(id='reaction_0'),
                'reaction_1': Reaction(id='reaction_1'),
                'reaction_2': Reaction(id='reaction_2'),
            },
        }

        value = "2*dfba_obj_reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']])
        self.assertEqual(invalid_attribute, None)
        rv = of_expr.validate()
        self.assertEqual(rv, None, str(rv))

        value = "2*dfba_obj_reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(dfba_obj_reactions=[])
        self.assertEqual(invalid_attribute, None)
        rv = of_expr.validate()
        self.assertNotEqual(rv, None, str(rv))

        value = "2*dfba_obj_reaction_1 - reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']],
            reactions=objs[Reaction].values())
        self.assertEqual(invalid_attribute, None)
        rv = of_expr.validate()
        self.assertEqual(rv, None, str(rv))

        value = "2*dfba_obj_reaction_1 - reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']],
            reactions=[])
        self.assertEqual(invalid_attribute, None)
        rv = of_expr.validate()
        self.assertNotEqual(rv, None, str(rv))

        value = "2*dfba_obj_reaction_1 - reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.expression = of_expr.expression[0:-1]
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']],
            reactions=objs[Reaction].values())
        rv = of_expr.validate()
        self.assertIsInstance(rv, InvalidObject)
        self.assertRegex(rv.attributes[0].messages[0], re.escape("aren't the id(s) of an object"))

        value = "2*dfba_obj_reaction_1 - reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.expression += ')'
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']],
            reactions=objs[Reaction].values())
        rv = of_expr.validate()
        self.assertIsInstance(rv, InvalidObject)
        self.assertRegex(rv.attributes[0].messages[0], "Python syntax error")

        value = "2*dfba_obj_reaction_1 -  3*reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.dfba_obj_reactions = []
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']],
            reactions=objs[Reaction].values())
        rv = of_expr.validate()
        self.assertRegex(rv.attributes[0].messages[0], re.escape("aren't the id(s) of an object"))

        value = "2*dfba_obj_reaction_1 * reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']],
            reactions=objs[Reaction].values())
        self.assertEqual(invalid_attribute, None)
        rv = of_expr.validate()
        self.assertEqual(rv, None)
        self.assertFalse(of_expr._parsed_expression.is_linear)

        value = "2*dfba_obj_reaction_1 ** 2"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']],
            reactions=objs[Reaction].values())
        self.assertEqual(invalid_attribute, None)
        rv = of_expr.validate()
        self.assertEqual(rv, None)
        self.assertFalse(of_expr._parsed_expression.is_linear)

        value = "2*dfba_obj_reaction_1 - pow( reaction_1, 2)"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        self.assertNotEqual(invalid_attribute, None)

        value = "1. + dfba_obj_reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']])
        rv = of_expr.validate()
        self.assertEqual(rv, None)
        self.assertFalse(of_expr._parsed_expression.is_linear)

        value = "1 + dfba_obj_reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']])
        rv = of_expr.validate()
        self.assertEqual(rv, None)
        self.assertFalse(of_expr._parsed_expression.is_linear)

        of_expr = DfbaObjectiveExpression(expression="True + dfba_obj_reaction_1",
                                          dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']])
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']])
        rv = of_expr.validate()
        self.assertRegex(rv.attributes[0].messages[0], re.escape("which aren't the id(s) of an object"))

        of_expr = DfbaObjectiveExpression(expression="'str' + dfba_obj_reaction_1",
                                          dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']])
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(
            dfba_obj_reactions=[objs[DfbaObjReaction]['dfba_obj_reaction_1']],
            reactions=objs[Reaction].values())
        rv = of_expr.validate()
        self.assertRegex(rv.attributes[0].messages[0], re.escape("contains bad token(s)"))

        of = DfbaObjective(id='dfba-obj-submdl_1',
                           submodel=Submodel(id='submdl_1'),
                           expression=DfbaObjectiveExpression(expression='1.'))
        self.assertEqual(of.validate(), None)

        of = DfbaObjective(id='invalid_id',
                           submodel=Submodel(id='submdl_2'),
                           expression=DfbaObjectiveExpression(expression='1.'))
        self.assertNotEqual(of.validate(), None)

        of = DfbaObjective(id='dfba-obj-submdl_3',
                           expression=DfbaObjectiveExpression(expression='1.'))
        self.assertNotEqual(of.validate(), None)

        value = "3*dfba_obj_reaction_1 - reaction_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value, objs)
        of = DfbaObjective(id='dfba-obj-submdl_4',
                           submodel=Submodel(id='submdl_4'),
                           expression=of_expr)
        of.submodel.dfba_obj_reactions = [objs[DfbaObjReaction]['dfba_obj_reaction_1']]
        of.submodel.reactions = objs[Reaction].values()
        errors = of_expr.validate()
        self.assertEqual(errors, None, str(errors))
        of_expr.reactions.append(objs[Reaction]['reaction_0'])
        self.assertNotEqual(of_expr.validate(), None)

        dfba_obj_rxn_1 = DfbaObjReaction(id='dfba_obj_rxn_1')
        value = "dfba_obj_rxn_1"
        of_expr, invalid_attribute = DfbaObjectiveExpression.deserialize(value,
                                                                         {DfbaObjReaction: {'dfba_obj_rxn_1': dfba_obj_rxn_1}})
        of_expr.dfba_obj = DfbaObjective()
        of_expr.dfba_obj.submodel = Submodel(dfba_obj_reactions=[dfba_obj_rxn_1])
        self.assertEqual(invalid_attribute, None)
        rv = of_expr.validate()
        self.assertEqual(rv, None, str(rv))

        of_expr.reactions = [dfba_obj_rxn_1]
        rv = of_expr.validate()
        self.assertNotEqual(rv, None, str(rv))

    def test_dfba_obj_species_validate(self):
        dfba_obj_reaction = DfbaObjReaction(id='dfba_obj_reaction', cell_size_units=unit_registry.parse_units('l'))
        species = Species(id='species')
        comp = DfbaObjSpecies(dfba_obj_reaction=dfba_obj_reaction,
                              species=species)
        comp.id = comp.gen_id()

        rv = comp.validate()
        self.assertEqual(rv, None, str(rv))

        comp.id = 'a' * 1000
        rv = comp.validate()
        self.assertNotEqual(rv, None, str(rv))

        comp.id = 'dfba-net-species-' + 'dfba_obj_reaction' + '_' + 'species'
        rv = comp.validate()
        self.assertNotEqual(rv, None, str(rv))

        comp.id = 'dfba-net-species-' + 'dfba_obj_reaction' + '-' + 'species'
        comp.units = unit_registry.parse_units('mol gDCW^-1 s^-1')
        rv = comp.validate()
        self.assertNotEqual(rv, None, str(rv))

    def test_validate(self):
        invalid = self.model.validate()
        self.assertEqual(invalid, None, str(invalid))

        self.model.id = 'invalid-id'
        invalid = self.model.validate()
        self.assertNotEqual(invalid, None, str(invalid))

    def test_dfba_obj_get_products(self):
        model = Model()
        submodel = model.submodels.create()
        species_type_0 = model.species_types.create(id='spec_0')
        species_type_1 = model.species_types.create(id='spec_1')
        species_type_2 = model.species_types.create(id='spec_2')
        species_type_3 = model.species_types.create(id='spec_3')
        compartment_0 = model.compartments.create(id='c_0')
        compartment_1 = model.compartments.create(id='c_1')
        compartment_2 = model.compartments.create(id='c_2')
        compartment_3 = model.compartments.create(id='c_3')
        species_0 = Species(id='spec_0[c_0]', species_type=species_type_0, compartment=compartment_0)
        species_1 = Species(id='spec_1[c_1]', species_type=species_type_1, compartment=compartment_1)
        species_2 = Species(id='spec_2[c_2]', species_type=species_type_2, compartment=compartment_2)
        species_3 = Species(id='spec_3[c_3]', species_type=species_type_3, compartment=compartment_3)

        obj_func = submodel.dfba_obj = DfbaObjective(
            expression=DfbaObjectiveExpression(
                reactions=[
                    Reaction(
                        reversible=True,
                        participants=[SpeciesCoefficient(species=species_0)],
                    ),
                    Reaction(
                        reversible=False,
                        participants=[
                            SpeciesCoefficient(species=species_1, coefficient=-1),
                            SpeciesCoefficient(species=species_2, coefficient=1),
                        ],
                    ),
                ],
                dfba_obj_reactions=[
                    DfbaObjReaction(
                        id='dfba_obj_rxn',
                        dfba_obj_species=[
                            DfbaObjSpecies(value=-1, species=species_1),
                            DfbaObjSpecies(value=1, species=species_3),
                        ],
                    ),
                ],
            ),
        )
        for rxn in obj_func.expression.dfba_obj_reactions:
            for spec in rxn.dfba_obj_species:
                spec.id = spec.gen_id()

        self.assertEqual(set(obj_func.get_products()), set([species_0, species_2, species_3]))
        self.assertEqual(set(obj_func.get_products(__type=Species)), set([species_0, species_2, species_3]))
        self.assertEqual(obj_func.get_products(__type=Model), [])

    def make_objects(self):
        model = Model()
        objects = {
            Observable: {},
            Parameter: {},
            Function: {},
            Species: {},
        }
        for id in ['a', 'b', 'duped_id']:
            param = model.parameters.create(id=id, value=1.)
            objects[Parameter][id] = param

        # use existing species
        for s in self.species:
            objects[Species][s.id] = s

        for id, species in zip(['ccc', 'ddd', 'eee', 'duped_id'], self.species):
            observable = model.observables.create(id=id)
            observable.expression, _ = ObservableExpression.deserialize(species.id, objects)
            objects[Observable][id] = observable

        for id, species in zip(['f', 'g', 'duped_id'], self.species):
            function = model.functions.create(id=id)
            function.expression, _ = FunctionExpression.deserialize(species.id, objects)
            objects[Function][id] = function

        id_map = {}
        for model_type in objects.keys():
            for id, obj in objects[model_type].items():
                typed_id = "{}.{}".format(model_type.__name__, id)
                id_map[typed_id] = obj

        return model, objects, id_map

    def test_make_obj(self):
        model, objects, id_map = self.make_objects()
        expr = 'ccc + 2 * ddd'
        fun_obj = Expression.make_obj(model, Function, 'fun_id', expr, objects)
        self.assertIsInstance(fun_obj, Function)
        self.assertEqual(fun_obj.id, 'fun_id')
        self.assertEqual(fun_obj.model, model)
        self.assertIsInstance(fun_obj.expression, FunctionExpression)
        self.assertTrue(fun_obj in model.functions)
        fun_expr = fun_obj.expression
        self.assertEqual(fun_expr.expression, expr)
        self.assertEqual(set(fun_expr.observables),
                         set([id_map['Observable.ccc'], id_map['Observable.ddd']]))
        self.assertEqual(set(fun_expr.parameters), set([]))
        self.assertEqual(set(fun_expr.functions), set([]))

        expr = 'ccc + 2 * x'
        expr_model_obj, error = Expression.make_expression_obj(Function, expr, objects)
        self.assertTrue(expr_model_obj is None)
        self.assertIsInstance(error, InvalidAttribute)

        fun_obj = Expression.make_obj(model, Function, 'fun_id', '', objects,
                                      allow_invalid_objects=True)
        self.assertNotIsInstance(fun_obj, Function)

    def test_make_obj_error(self):
        model, objects, id_map = self.make_objects()

        obs = Expression.make_obj(model, Observable, 'obs_id', 'ccc + ddd', objects)
        self.assertIsInstance(obs, Observable)

        error = Expression.make_obj(model, Observable, 'obs_id', 'ccc * ddd', objects)
        self.assertNotIsInstance(error, Observable)
        self.assertIsInstance(error, InvalidObject)

    def do_test_valid_expression(self, expression_class, parent_class, objects, expr, expected_val,
                                 expected_related_objs=None, expected_error=None):
        """ Test a valid expression

        Args:
            expression_class (:obj:`obj_tables.Model`): expression class being tested
            parent_class (:obj:`obj_tables.Model`): the expression model that uses an expression_class
            objects (:obj:`dict`): dict of objects which can be used by `expr`
            expr (:obj:`str`): the expression
            expected_val (:obj:`obj`): the value expected when the expression is evaluated
            expected_related_objs (:obj:`dict`, optional): objects that should be used by the deserialize expression
        """
        expr_obj, error = expression_class.deserialize(expr, objects)
        parent = parent_class(expression=expr_obj)
        if expected_error:
            self.assertIn(expected_error, str(error))
        else:
            self.assertEqual(error, None, str(error))
            self.assertEqual(expr_obj.serialize(), expr)
            self.assertEqual(expr_obj.expression, expr)
            self.assertEqual(expr_obj, objects[expression_class][expr])
            # check used objects in expression_class attributes
            if expected_related_objs:
                for modifier, elements in expected_related_objs.items():
                    self.assertEqual(set(getattr(expr_obj, modifier)), set(elements))
            error = expr_obj.validate()
            self.assertEqual(error, None, str(error))
            self.assertEqual(expr_obj._parsed_expression.test_eval({Species: 1.}),
                             expected_val,
                             expr_obj._parsed_expression.expression)

    def test_valid_function_expressions(self):
        _, objects, id_map = self.make_objects()

        i = 0
        for i_expr, (expr, expected_test_val, expected_related_objs, error) in enumerate([
            ('ccc', 1, {'observables': [id_map['Observable.ccc']]}, None),
            ('ddd + eee', 2, {'observables': [id_map['Observable.ddd'], id_map['Observable.eee']]}, None),
            ('ddd + 2 * eee', 3, {}, None),
            ('ddd + 2 * eee > 3', False, {}, None),
            ('1 * duped_id', None, {}, 'multiple model object id matches'),
            ('1 * Observable.duped_id', 1,
                {'observables': [id_map['Observable.duped_id']]},
                None),
            ('a + f', 2,
                {'functions': [id_map['Function.f']]},
                None),
            ('log(a)', math.log(1), {}, None),
            ('max(a, b)', 1,
                {'parameters': [id_map['Parameter.a'], id_map['Parameter.b']]},
                None),
            ('Observable.ddd + Observable.duped_id - Parameter.duped_id', 1,
                {'parameters': [id_map['Parameter.duped_id']],
                 'observables':[id_map['Observable.duped_id'], id_map['Observable.ddd']]},
                None),
            ('ddd * Function.duped_id', 1,
                {'observables': [id_map['Observable.ddd']],
                 'functions':[id_map['Function.duped_id']]},
                None),
        ]):
            self.do_test_valid_expression(FunctionExpression, Function,
                                          objects, expr, expected_test_val, expected_related_objs,
                                          expected_error=error)

        # reuse the FunctionExpression
        with self.assertRaisesRegex(ValueError, 'must be `None`'):
            self.do_test_valid_expression(FunctionExpression, Function,
                                          objects, 'ccc', 1, {'observables': [id_map['Observable.ccc']]},
                                          expected_error=None)

    def do_test_expr_deserialize_error(self, expression_class, parent_class, objects, expr, error_msg_substr):
        """ Test an expression that fails to deserialize

        Args:
            expression_class (:obj:`obj_tables.Model`): expression class being tested
            parent_class (:obj:`obj_tables.Model`): the expression model that uses an expression_class
            objects (:obj:`dict`): dict of objects which can be used by `expr`
            expr (:obj:`str`): the expression
            error_msg_substr (:obj:`str`): substring expected in error message
        """
        func_expr, error = expression_class.deserialize(expr, objects)
        self.assertEqual(func_expr, None)
        self.assertIsInstance(error, InvalidAttribute)
        self.assertIn(error_msg_substr, error.messages[0])

    def test_function_expression_deserialize_errors(self):
        _, objects, _ = self.make_objects()

        for expr_model, model in [(FunctionExpression, Function),
                                  (StopConditionExpression, StopCondition),
                                  (ObservableExpression, Observable)]:
            self.do_test_expr_deserialize_error(expr_model, model, objects, 'id1[id2', "Python syntax error")
            bad_id = 'no_such_obj'
            self.do_test_expr_deserialize_error(expr_model, model, objects, bad_id,
                                                "contains the identifier(s) '{}', which aren't the id(s) of an object".format(bad_id))

    def test_stop_condition_expression_deserialize_errors(self):
        _, objects, _ = self.make_objects()

        expr = '(ccc > 10 and ddd < 5) or (a + f * g)'
        self.do_test_expr_deserialize_error(StopConditionExpression, StopCondition, objects, expr,
                                            "contains the identifier(s) 'and', which aren't the id(s) of an object")

    def do_test_invalid_expression(self, expression_class, parent_class, objects, expr, error_msg_substr):
        """ Test an expression that fails to validate

        Args:
            expression_class (:obj:`obj_tables.Model`): expression class being tested
            parent_class (:obj:`obj_tables.Model`): the expression model that uses an expression_class
            objects (:obj:`dict`): dict of objects which can be used by `expr`
            expr (:obj:`str`): the expression
            error_msg_substr (:obj:`str`): substring expected in error message
        """
        func_expr, error = expression_class.deserialize(expr, objects)
        self.assertEqual(error, None)
        invalid_obj = func_expr.validate()
        self.assertIsInstance(invalid_obj, InvalidObject)
        self.assertIn(error_msg_substr, invalid_obj.attributes[0].messages[0])

    def test_invalid_function_expressions(self):
        func_expr, error = FunctionExpression.deserialize('1 +', {})
        self.assertRegex(str(error), "SyntaxError:")

    def test_function(self):
        model, objects, _ = self.make_objects()
        kwargs = dict(id='fun_1', name='name fun_1', model=model, comments='no comment')
        func = Function(**kwargs)
        self.assertEqual(func.id, kwargs['id'])
        self.assertEqual(func.name, kwargs['name'])
        self.assertEqual(func.model, model)
        self.assertEqual(model.functions[-1], func)
        self.assertEqual(func.comments, kwargs['comments'])
        expr = 'ccc + ddd'
        func_expr, _ = FunctionExpression.deserialize(expr, objects)
        func.expression = func_expr

    def test_function_validate(self):
        func = Function(id='func_0', units=unit_registry.parse_units('dimensionless'))
        func_expr, _ = FunctionExpression.deserialize('1', {})
        func.expression = func_expr
        self.assertEqual(func.validate(), None)

        func = Function(id='func-0', units=unit_registry.parse_units('dimensionless'))
        func_expr, _ = FunctionExpression.deserialize('1', {})
        func.expression = func_expr
        self.assertNotEqual(func.validate(), None)

    def test_function_deserialize_invalid_ids(self):

        objs = {
            SpeciesType: {
                'st_1': SpeciesType(id='st_1'),
            },
            Compartment: {
                'c': Compartment(id='c'),
            },
            Species: {},
            Parameter: {
                'pow': Parameter(id='pow'),
            }
        }
        objs[Species]['st_1[c]'] = Species(id='st_1[c]',
                                           species_type=objs[SpeciesType]['st_1'],
                                           compartment=objs[Compartment]['c'])

        value = "pow( st_1[c], 2 )"
        of_expr, invalid_attribute = FunctionExpression.deserialize(value, objs)
        self.assertEqual(of_expr, None, str())
        self.assertIn("ObjTablesToken `pow` is ambiguous",
                      invalid_attribute.messages[0])
        self.assertIn("ObjTablesToken matches a Parameter and a function",
                      invalid_attribute.messages[0])

    def test_valid_stop_conditions(self):
        _, objects, id_map = self.make_objects()

        some_used_objs = {'observables': [id_map['Observable.ccc'], id_map['Observable.ddd']],
                          'functions': [id_map['Function.f'], id_map['Function.g']],
                          'parameters': [id_map['Parameter.a']]}
        for expr, expected_test_val, expected_attrs in [
            ('ccc > 10', False, {'observables': [id_map['Observable.ccc']]}),
            ('ccc > 0', True, {'observables': [id_map['Observable.ccc']]}),
            ('ccc + ddd - a + f * g + 10 > 0', True, some_used_objs)
        ]:
            self.do_test_valid_expression(StopConditionExpression, StopCondition,
                                          objects, expr, expected_test_val, expected_attrs)

        with self.assertRaisesRegex(ValueError, 'must be `None`'):
            # reuse StopConditionExpression
            self.do_test_valid_expression(StopConditionExpression, StopCondition,
                                          objects, 'ccc > 0', True, {'observables': [id_map['Observable.ccc']]})

    def test_invalid_stop_condition_expressions(self):
        _, objects, _ = self.make_objects()

        bad_expr = '1 + ccc'
        self.do_test_invalid_expression(StopConditionExpression, StopCondition, objects, bad_expr,
                                        "Evaluating '{}', a {} expression, should return a bool but it returns a float".format(
                                            bad_expr, StopConditionExpression.__name__))

    def test_stop_condition(self):
        model, objects, _ = self.make_objects()
        kwargs = dict(id='stop_cond_1', name='name stop_cond_1', model=model, comments='no comment')
        stop_condition = StopCondition(**kwargs)
        self.assertEqual(stop_condition.id, kwargs['id'])
        self.assertEqual(stop_condition.name, kwargs['name'])
        self.assertEqual(stop_condition.model, model)
        self.assertEqual(model.stop_conditions[-1], stop_condition)
        self.assertEqual(stop_condition.comments, kwargs['comments'])

        expr = 'ccc + ddd'
        stop_condition_expr, _ = StopConditionExpression.deserialize(expr, objects)
        stop_condition.expression = stop_condition_expr
        self.assertNotEqual(stop_condition.validate(), None)

        expr = 'ccc / ddd > 2'
        stop_condition_expr, _ = StopConditionExpression.deserialize(expr, objects)
        stop_condition.expression = stop_condition_expr
        rv = stop_condition.validate()
        self.assertEqual(rv, None, str(rv))

        expr = 'ccc / !ddd > 2'
        stop_condition_expr, _ = StopConditionExpression.deserialize(expr, objects)
        stop_condition.expression = stop_condition_expr
        self.assertEqual(stop_condition_expr, None)
        rv = stop_condition.validate()
        self.assertNotEqual(rv, None, str(rv))

        expr = 'ccc / ~ddd > 2'
        stop_condition_expr, _ = StopConditionExpression.deserialize(expr, objects)
        stop_condition.expression = stop_condition_expr
        self.assertEqual(stop_condition_expr, None)
        rv = stop_condition.validate()
        self.assertNotEqual(rv, None, str(rv))

        expr = 'ccc / ddd > 2\n'
        stop_condition_expr, _ = StopConditionExpression.deserialize(expr, objects)
        self.assertNotEqual(stop_condition_expr, None)

        expr = 'ccc << ddd'
        stop_condition_expr, _ = StopConditionExpression.deserialize(expr, objects)
        self.assertEqual(stop_condition_expr, None)

    def test_stop_condition_validate(self):
        stop_cond = StopCondition(id='stop_cond_0', units=unit_registry.parse_units('dimensionless'))
        stop_cond_expr, _ = StopConditionExpression.deserialize('2 > 1', {})
        stop_cond.expression = stop_cond_expr
        self.assertEqual(stop_cond.validate(), None)

        stop_cond = StopCondition(id='stop_cond-0', units=unit_registry.parse_units('dimensionless'))
        stop_cond_expr, _ = StopConditionExpression.deserialize('2 > 1', {})
        stop_cond.expression = stop_cond_expr
        self.assertNotEqual(stop_cond.validate(), None)

        stop_cond = StopCondition(id='stop_cond_0', units=unit_registry.parse_units('s'))
        stop_cond_expr, _ = StopConditionExpression.deserialize('2 > 1', {})
        stop_cond.expression = stop_cond_expr
        self.assertNotEqual(stop_cond.validate(), None)

        stop_cond = StopCondition(id='stop_cond_0', units=unit_registry.parse_units('dimensionless'))
        stop_cond_expr, _ = StopConditionExpression.deserialize('2 > 1', {})
        stop_cond.expression = stop_cond_expr
        stop_cond.expression._parsed_expression = None
        self.assertNotEqual(stop_cond.validate(), None)

    def test_valid_observable_expressions(self):
        _, objects, id_map = self.make_objects()

        for expr, expected_test_val, expected_related_objs in [
            ('3 * spec_type_0[comp_0]', 3, {'species': [id_map['Species.spec_type_0[comp_0]']]}),
            ('spec_type_0[comp_0] + 4e2*spec_type_1[comp_0] - 1  * spec_type_3[comp_1]', 400,
                {'species': [
                    id_map['Species.spec_type_0[comp_0]'],
                    id_map['Species.spec_type_1[comp_0]'],
                    id_map['Species.spec_type_3[comp_1]']
                ]}),
            ('1.5 * spec_type_0[comp_0] - ddd + Observable.duped_id', 1.5,
                {'species': [id_map['Species.spec_type_0[comp_0]']],
                 'observables':[id_map['Observable.duped_id'], id_map['Observable.ddd']],
                 }),
        ]:
            self.do_test_valid_expression(ObservableExpression, Observable,
                                          objects, expr, expected_test_val, expected_related_objs)

    def test_invalid_observable_expressions(self):
        _, objects, _ = self.make_objects()

        non_linear_expression = 'ccc * ccc'
        self.do_test_invalid_expression(ObservableExpression, Observable, objects, non_linear_expression,
                                        "Expression must be linear")

    def test_observable(self):
        model, objects, _ = self.make_objects()
        kwargs = dict(id='obs_1', name='name obs_1', model=model, comments='no comment')
        obs = Observable(**kwargs)
        self.assertEqual(obs.id, kwargs['id'])
        self.assertEqual(obs.name, kwargs['name'])
        self.assertEqual(obs.model, model)
        self.assertEqual(model.observables[-1], obs)
        self.assertEqual(obs.comments, kwargs['comments'])

        expr = 'ccc + ddd - 2 * spec_type_0[comp_0]'
        obs_expr, _ = ObservableExpression.deserialize(expr, objects)
        self.assertEqual(obs_expr.validate(), None)

        expr = 'ccc + ddd + 2 - 2 * spec_type_0[comp_0]'
        obs_expr, _ = ObservableExpression.deserialize(expr, objects)
        self.assertNotEqual(obs_expr.validate(), None)

        expr = 'ccc + ddd - spec_type_0[comp_0] * 2'
        obs_expr, _ = ObservableExpression.deserialize(expr, objects)
        self.assertNotEqual(obs_expr.validate(), None)

        expr = 'ccc'
        obs_expr, _ = ObservableExpression.deserialize(expr, objects)
        self.assertEqual(obs_expr.validate(), None)

        expr = 'ccc * ccc'
        obs_expr, _ = ObservableExpression.deserialize(expr, objects)
        self.assertNotEqual(obs_expr.validate(), None)

        expr = '2'
        obs_expr, _ = ObservableExpression.deserialize(expr, objects)
        self.assertNotEqual(obs_expr.validate(), None)

        self.assertEqual(Observable.Meta.attributes['expression'].serialize(None), '')
        self.assertEqual(Observable.Meta.attributes['expression'].serialize(''), '')
        self.assertEqual(Observable.Meta.attributes['expression'].deserialize(None, objects), (None, None))
        self.assertEqual(Observable.Meta.attributes['expression'].deserialize('', objects), (None, None))

    def test_expression_term_model_types(self):
        for model_type in [RateLawExpression, FunctionExpression, StopConditionExpression,
                           DfbaObjectiveExpression, ObservableExpression]:
            self.assertTrue(hasattr(model_type.Meta, 'expression_term_models'))
            for expression_term_model_type in model_type.Meta.expression_term_models:
                self.assertTrue(hasattr(wc_lang.core, expression_term_model_type))

    def test_reaction_get_reactants(self):
        rxn = Reaction()
        species_1 = Species()
        species_2 = Species()
        species_3 = Species()
        species_4 = Species()
        species_5 = Species()
        rxn.participants.append(SpeciesCoefficient(coefficient=1., species=species_1))
        rxn.participants.append(SpeciesCoefficient(coefficient=2, species=species_2))
        rxn.participants.append(SpeciesCoefficient(coefficient=-1., species=species_3))
        rxn.participants.append(SpeciesCoefficient(coefficient=0, species=species_4))
        rxn.participants.append(SpeciesCoefficient(coefficient=-2, species=species_5))

        self.assertEqual(rxn.get_reactants(), [species_3, species_5])

    def test_reaction_get_products(self):
        rxn = Reaction()
        species_1 = Species()
        species_2 = Species()
        species_3 = Species()
        species_4 = Species()
        species_5 = Species()
        rxn.participants.append(SpeciesCoefficient(coefficient=1., species=species_1))
        rxn.participants.append(SpeciesCoefficient(coefficient=2, species=species_2))
        rxn.participants.append(SpeciesCoefficient(coefficient=-1., species=species_3))
        rxn.participants.append(SpeciesCoefficient(coefficient=0, species=species_4))
        rxn.participants.append(SpeciesCoefficient(coefficient=-2, species=species_5))

        self.assertEqual(rxn.get_products(), [species_1, species_2])

    def test_reaction_get_modifiers(self):
        rxn = Reaction()
        species_1 = Species()
        species_2 = Species()
        species_3 = Species()
        species_4 = Species()
        species_5 = Species()
        species_6 = Species()
        rxn.participants.append(SpeciesCoefficient(coefficient=1., species=species_1))
        rxn.participants.append(SpeciesCoefficient(coefficient=2, species=species_2))
        rxn.participants.append(SpeciesCoefficient(coefficient=-1., species=species_3))
        rxn.participants.append(SpeciesCoefficient(coefficient=0, species=species_4))
        rxn.participants.append(SpeciesCoefficient(coefficient=-2, species=species_5))

        rl = rxn.rate_laws.create()
        rl.expression = RateLawExpression()
        rl.expression.species = [species_3, species_6, species_1, species_5]
        self.assertEqual(set(rxn.get_modifiers()), set([species_6, species_1]))

    def test_author_get_identifier(self):
        author = Author()
        author.identifiers.create(namespace='github.user', id='jonrkarr')
        author.identifiers.create(namespace='orcid', id='0000-0002-2605-5080')
        author.identifiers.create(namespace='github.organization', id='karrlab')
        author.identifiers.create(namespace='github.organization', id='wholecell')
        self.assertEqual(author.get_identifier('github.user'), 'jonrkarr')
        self.assertEqual(author.get_identifier('orcid'), '0000-0002-2605-5080')
        self.assertEqual(author.get_identifier('linkedin.user'), None)
        with self.assertRaisesRegex(ValueError, 'has multiple'):
            author.get_identifier('github.organization')

    def test_ChemicalStructure_get_structure(self):
        s = ChemicalStructure()
        self.assertEqual(s.get_structure(), None)

        s = ChemicalStructure(value='[OH2]', format=ChemicalStructureFormat.SMILES)
        self.assertEqual(s.get_structure().GetTotalCharge(), 0)

        s = ChemicalStructure(value='AAA', format=ChemicalStructureFormat.BpForms, alphabet=ChemicalStructureAlphabet.dna)
        self.assertEqual(s.get_structure().get_charge(), -4)

        s = ChemicalStructure(value='AAA', format='BpForms')
        with self.assertRaisesRegex(ValueError, 'Unsupported format'):
            s.get_structure()

    def test_ChemicalStructure_validate(self):
        s = ChemicalStructure()
        err = s.validate()
        self.assertEqual(err, None, str(err))

        s = ChemicalStructure(molecular_weight=1., charge=1)
        err = s.validate()
        self.assertEqual(err, None, str(err))

        s = ChemicalStructure(molecular_weight=1.)
        err = s.validate()
        self.assertEqual(err, None, str(err))

        s = ChemicalStructure(charge=1)
        err = s.validate()
        self.assertEqual(err, None, str(err))

        s = ChemicalStructure(empirical_formula=EmpiricalFormula('OH'), charge=-1)
        err = s.validate()
        self.assertEqual(err, None, str(err))

        s = ChemicalStructure(empirical_formula=EmpiricalFormula('OH'), charge=-1)
        s.molecular_weight = s.empirical_formula.get_molecular_weight()
        err = s.validate()
        self.assertEqual(err, None, str(err))

        s = ChemicalStructure(value='[OH2]', format=ChemicalStructureFormat.SMILES,
                              empirical_formula=EmpiricalFormula('OH2'), charge=0)
        s.molecular_weight = s.empirical_formula.get_molecular_weight()
        err = s.validate()
        self.assertEqual(err, None, str(err))

        s = ChemicalStructure(value='AAA', format=ChemicalStructureFormat.BpForms,
                              alphabet=ChemicalStructureAlphabet.dna)
        form = s.get_structure()
        s.empirical_formula = form.get_formula()
        s.molecular_weight = form.get_mol_wt()
        s.charge = form.get_charge()
        err = s.validate()
        self.assertEqual(err, None, str(err))

        s = ChemicalStructure(value='AAA', alphabet=ChemicalStructureAlphabet.dna)
        err = s.validate()
        self.assertNotEqual(err, None, str(err))

        s = ChemicalStructure(format=ChemicalStructureFormat.BpForms,
                              alphabet=ChemicalStructureAlphabet.dna)
        err = s.validate()
        self.assertNotEqual(err, None, str(err))

        s = ChemicalStructure(value='AAA', format=ChemicalStructureFormat.BpForms)
        err = s.validate()
        self.assertNotEqual(err, None, str(err))

        s = ChemicalStructure(value='[OH2]', format=ChemicalStructureFormat.SMILES,
                              alphabet=ChemicalStructureAlphabet.dna)
        err = s.validate()
        self.assertNotEqual(err, None, str(err))

        s = ChemicalStructure(value='[OH2]', format=ChemicalStructureFormat.SMILES,
                              empirical_formula=EmpiricalFormula('H2O'), charge=0)
        mol_wt = s.molecular_weight = s.empirical_formula.get_molecular_weight()
        err = s.validate()
        self.assertEqual(err, None, str(err))

        s = ChemicalStructure(value='[OH2]', format=ChemicalStructureFormat.SMILES,
                              empirical_formula=EmpiricalFormula('H3O'), charge=0, molecular_weight=mol_wt)
        err = s.validate()
        self.assertNotEqual(err, None, str(err))

        s = ChemicalStructure(value='[OH2]', format=ChemicalStructureFormat.SMILES,
                              empirical_formula=EmpiricalFormula('H2O'), charge=-1, molecular_weight=mol_wt)
        err = s.validate()
        self.assertNotEqual(err, None, str(err))

        s = ChemicalStructure(value='[OH2]', format=ChemicalStructureFormat.SMILES,
                              empirical_formula=EmpiricalFormula('H2O'), charge=0, molecular_weight=mol_wt + 1)
        err = s.validate()
        self.assertNotEqual(err, None, str(err))

        s = ChemicalStructure(value='[OH2]', format=ChemicalStructureFormat.SMILES,
                              empirical_formula=EmpiricalFormula('H2O'), charge=0, molecular_weight=-1.)
        err = s.validate()
        self.assertNotEqual(err, None, str(err))

    def test_evidence_serialize_deserialize(self):
        obs_1 = Observation(id='obs_1')
        ev = Evidence(observation=obs_1, type=onto['WC:supporting_evidence'], strength=10., quality=20.)
        self.assertEqual(ev.serialize(), 'obs_1(+, s=10.0, q=20.0)')

        obs_1b = Observation(id='obs_1b')
        ev = Evidence(observation=obs_1b, type=onto['WC:supporting_evidence'], strength=10., quality=20.)
        self.assertEqual(ev.serialize(), 'obs_1b(+, s=10.0, q=20.0)')

        obs_2 = Observation(id='obs_2')
        ev = Evidence(observation=obs_2, type=onto['WC:disputing_evidence'], strength=10.)
        self.assertEqual(ev.serialize(), 'obs_2(-, s=10.0)')

        obs_3 = Observation(id='obs_3')
        ev = Evidence(observation=obs_3, type=onto['WC:inconclusive_evidence'], quality=20.)
        self.assertEqual(ev.serialize(), 'obs_3(~, q=20.0)')

        obs_4 = Observation(id='obs_4')
        ev = Evidence(observation=obs_4, type=onto['WC:supporting_evidence'])
        self.assertEqual(ev.serialize(), 'obs_4(+)')

        attr = Submodel.Meta.attributes['evidence']
        obs_7 = Observation(id='obs_7')
        ev = [
            Evidence(observation=obs_7, type=onto['WC:supporting_evidence'], strength=10.),
            Evidence(observation=obs_7, type=onto['WC:supporting_evidence'], strength=20.),
        ]
        self.assertEqual(attr.serialize(ev), 'obs_7(+, s=10.0); obs_7(+, s=20.0)')


class TestCoreFromFile(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_model.xlsx')

    def setUp(self):
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)[Model][0]
        self.dfba_submodel = self.model.submodels.get_one(id='submodel_1')

    def test_get_species(self):
        species_ids = set([s.id for s in self.dfba_submodel.get_children(kind='submodel', __type=Species)])
        self.assertEqual(species_ids, set([
            'specie_1[c]',
            'specie_1[e]',
            'specie_2[c]',
            'specie_2[e]',
            'specie_3[c]',
            'specie_4[c]',
        ]))


class TestTaxonRank(unittest.TestCase):

    def test(self):
        self.assertEqual(TaxonRank['varietas'], TaxonRank.variety)
        self.assertEqual(TaxonRank['strain'], TaxonRank.variety)
        self.assertEqual(TaxonRank['tribus'], TaxonRank.tribe)
        self.assertEqual(TaxonRank['familia'], TaxonRank.family)
        self.assertEqual(TaxonRank['ordo'], TaxonRank.order)
        self.assertEqual(TaxonRank['class'], TaxonRank.classis)
        self.assertEqual(TaxonRank['division'], TaxonRank.phylum)
        self.assertEqual(TaxonRank['divisio'], TaxonRank.phylum)
        self.assertEqual(TaxonRank['regnum'], TaxonRank.kingdom)


class ValidateModelTestCase(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_validate_model.xlsx')

    def setUp(self):
        # read a wc model
        self.model = Reader().run(self.MODEL_FILENAME)[Model][0]
        self.dfba_submodel = self.model.submodels.get_one(id='dfba_submodel')

    def test_min_flux_maxes(self):
        c_1 = Compartment(id='c_1')
        c_2 = Compartment(id='c_2')
        st = SpeciesType(id='s', structure=ChemicalStructure(empirical_formula=EmpiricalFormula('CHN2P1'), charge=-1))
        species_1 = Species(species_type=st, compartment=c_1)
        species_2 = Species(species_type=st, compartment=c_2)
        species_1.id = species_1.gen_id()
        species_2.id = species_2.gen_id()
        participants = [
            SpeciesCoefficient(species=species_1, coefficient=1.),
            SpeciesCoefficient(species=species_2, coefficient=-1.),
        ]

        rxn = Reaction(id='rxn', reversible=True,
                       flux_bounds=FluxBounds(min=-1., max=1., units=unit_registry.parse_units('M s^-1')),
                       submodel=Submodel(framework=onto['WC:dynamic_flux_balance_analysis']),
                       participants=participants)
        rv = rxn.validate()
        self.assertEqual(rv, None, str(rv))

        rxn = Reaction(id='rxn', reversible=False,
                       flux_bounds=FluxBounds(min=0., max=1., units=unit_registry.parse_units('M s^-1')),
                       submodel=Submodel(framework=onto['WC:dynamic_flux_balance_analysis']),
                       participants=participants)
        rv = rxn.validate()
        self.assertEqual(rv, None, str(rv))

        rxn = Reaction(id='rxn', reversible=True,
                       flux_bounds=FluxBounds(min=1., max=-1., units=unit_registry.parse_units('M s^-1')),
                       submodel=Submodel(framework=onto['WC:stochastic_simulation_algorithm']),
                       participants=participants,
                       rate_laws=[
                           RateLaw(direction=RateLawDirection.forward),
                           RateLaw(direction=RateLawDirection.backward),
                       ])
        rv = rxn.validate()
        self.assertEqual(len(rv.attributes), 3)
        self.assertRegex(str(rv), 'Flux bounds should be None')
        self.assertRegex(str(rv), 'Maximum flux must be least the minimum flux')
        self.assertRegex(str(rv), 'Minimum flux for reversible reaction should be negative')

        rxn = Reaction(id='rxn', reversible=False,
                       flux_bounds=FluxBounds(min=-1., max=-1.5, units=unit_registry.parse_units('M s^-1')),
                       submodel=Submodel(framework=onto['WC:stochastic_simulation_algorithm']),
                       participants=participants,
                       rate_laws=[
                           RateLaw(direction=RateLawDirection.forward),
                       ])
        rv = rxn.validate()
        self.assertEqual(len(rv.attributes), 3)
        self.assertRegex(str(rv), 'Flux bounds should be None')
        self.assertRegex(str(rv), 'Maximum flux must be least the minimum flux')
        self.assertRegex(str(rv), 'Minimum flux for irreversible reaction should be non-negative')

    def test_dfba_submodel_contains_obj_reactions(self):
        submodel = Submodel(id='submodel', framework=onto['WC:dynamic_flux_balance_analysis'])
        objs = {
            Reaction: {'rxn_1': Reaction(id='rxn_1', submodel=submodel)},
            DfbaObjReaction: {'dfba_obj_rxn_1': DfbaObjReaction(id='dfba_obj_rxn_1', submodel=submodel)}
        }
        of_expr, _ = DfbaObjectiveExpression.deserialize('rxn_1 + dfba_obj_rxn_1', objs)
        submodel.dfba_obj = of_expr.dfba_obj = DfbaObjective(id='submodel-dfba-obj')
        rv = submodel.validate()
        self.assertEqual(rv, None, str(rv))
        rv = of_expr.validate()
        self.assertEqual(rv, None, str(rv))

        submodel = Submodel(id='submodel', framework=onto['WC:dynamic_flux_balance_analysis'])
        rv = submodel.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertRegex(str(rv), 'submodel must have an objective')

        submodel = Submodel(id='submodel', framework=onto['WC:dynamic_flux_balance_analysis'])
        objs = {
            Reaction: {'rxn_1': Reaction(id='rxn_1', submodel=submodel)},
            DfbaObjReaction: {'dfba_obj_rxn_1': DfbaObjReaction(id='dfba_obj_rxn_1', submodel=submodel)}
        }
        of_expr, _ = DfbaObjectiveExpression.deserialize('rxn_1 + dfba_obj_rxn_1', objs)
        submodel.dfba_obj = of_expr.dfba_obj = DfbaObjective(id='submodel-dfba-obj')
        of_expr.reactions = []
        of_expr.dfba_obj_reactions = []
        rv = submodel.validate()
        self.assertEqual(rv, None, str(rv))
        rv = of_expr.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertEqual(len(rv.attributes[0].messages), 1)
        self.assertRegex(str(rv), 'must be a function of at least one')

        submodel = Submodel(id='submodel', framework=onto['WC:dynamic_flux_balance_analysis'])
        objs = {
            Reaction: {'rxn_1': Reaction(id='rxn_1', submodel=submodel)},
            DfbaObjReaction: {'dfba_obj_rxn_1': DfbaObjReaction(id='dfba_obj_rxn_1', submodel=submodel)}
        }
        of_expr, _ = DfbaObjectiveExpression.deserialize('rxn_1 + dfba_obj_rxn_1', objs)
        submodel.dfba_obj = of_expr.dfba_obj = DfbaObjective(id='submodel-dfba-obj')
        submodel.reactions = []
        submodel.dfba_obj_reactions = []
        rv = submodel.validate()
        self.assertEqual(rv, None, str(rv))
        rv = of_expr.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertEqual(len(rv.attributes[0].messages), 2)
        self.assertRegex(str(rv), 'must contain the following reactions')
        self.assertRegex(str(rv), 'must contain the following dFBA objective reactions')

        submodel = Submodel(id='submodel', framework=onto['WC:dynamic_flux_balance_analysis'])
        objs = {
            Reaction: {'rxn_1': Reaction(id='rxn_1', submodel=submodel)},
            DfbaObjReaction: {'dfba_obj_rxn_1': DfbaObjReaction(id='dfba_obj_rxn_1', submodel=submodel)}
        }
        of_expr, _ = DfbaObjectiveExpression.deserialize('rxn_1 + dfba_obj_rxn_1', objs)
        submodel.dfba_obj = of_expr.dfba_obj = DfbaObjective(id='submodel-dfba-obj')
        submodel.reactions = []
        rv = submodel.validate()
        self.assertEqual(rv, None, str(rv))
        rv = of_expr.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertEqual(len(rv.attributes[0].messages), 1)
        self.assertRegex(str(rv), 'must contain the following reactions')

    def test_rate_laws(self):
        c_1 = Compartment(id='c_1')
        c_2 = Compartment(id='c_2')
        st = SpeciesType(id='s', structure=ChemicalStructure(empirical_formula=EmpiricalFormula('CHO'), charge=1))
        species_1 = Species(species_type=st, compartment=c_1)
        species_2 = Species(species_type=st, compartment=c_2)
        species_1.id = species_1.gen_id()
        species_2.id = species_2.gen_id()
        participants = [
            SpeciesCoefficient(species=species_1, coefficient=-1.),
            SpeciesCoefficient(species=species_2, coefficient=1.),
        ]

        rxn = Reaction(id='rxn', reversible=True,
                       submodel=Submodel(framework=onto['WC:stochastic_simulation_algorithm']),
                       participants=participants,
                       rate_laws=[
                           RateLaw(direction=RateLawDirection.forward),
                           RateLaw(direction=RateLawDirection.backward),
                       ])
        rv = rxn.validate()
        self.assertEqual(rv, None, str(rv))

        rxn = Reaction(id='rxn', reversible=True,
                       submodel=Submodel(framework=onto['WC:stochastic_simulation_algorithm']),
                       participants=participants,
                       rate_laws=[
                       ])
        rv = rxn.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertEqual(len(rv.attributes[0].messages), 2)
        self.assertRegex(str(rv), 'must have a forward rate law')
        self.assertRegex(str(rv), 'must have a backward rate law')

        rxn = Reaction(id='rxn', reversible=False,
                       submodel=Submodel(framework=onto['WC:stochastic_simulation_algorithm']),
                       participants=participants,
                       rate_laws=[
                       ])
        rv = rxn.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertEqual(len(rv.attributes[0].messages), 1)
        self.assertRegex(str(rv), 'must have a forward rate law')

        rxn = Reaction(id='rxn', reversible=False,
                       submodel=Submodel(framework=onto['WC:stochastic_simulation_algorithm']),
                       participants=participants,
                       rate_laws=[
                           RateLaw(direction=RateLawDirection.forward),
                           RateLaw(direction=RateLawDirection.backward),
                       ])
        rv = rxn.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertEqual(len(rv.attributes[0].messages), 1)
        self.assertRegex(str(rv), 'cannot have a backward rate law')

    def test_species_types(self):
        st = SpeciesType(id='species_4', structure=ChemicalStructure(molecular_weight=1.))
        self.assertEqual(st.structure.validate(), None)

        st = SpeciesType(id='species_4', structure=ChemicalStructure(molecular_weight=0.))
        self.assertEqual(st.structure.validate(), None)

        st = SpeciesType(id='species_4', structure=ChemicalStructure(molecular_weight=-1.))
        self.assertNotEqual(st.structure.validate(), None)

    def test_acyclic_dependencies(self):
        model = Model(id='model', version='0.0.1')
        c = model.compartments.create(id='c')
        st_1 = model.species_types.create(id='st_1')
        st_2 = model.species_types.create(id='st_2')
        s_1 = model.species.create(species_type=st_1, compartment=c)
        s_2 = model.species.create(species_type=st_2, compartment=c)
        obs_1 = model.observables.create(id='obs_1',
                                         expression=ObservableExpression(expression='st_1[c]', species=[s_1]))
        obs_2 = model.observables.create(id='obs_2',
                                         expression=ObservableExpression(expression='st_2[c]', species=[s_2]))
        rv = model.validate()
        self.assertEqual(rv, None, str(rv))

        model = Model(id='model', version='0.0.1')
        c = model.compartments.create(id='c')
        st_1 = model.species_types.create(id='st_1')
        st_2 = model.species_types.create(id='st_2')
        s_1 = model.species.create(species_type=st_1, compartment=c)
        s_2 = model.species.create(species_type=st_2, compartment=c)
        obs_1 = model.observables.create(id='obs_1',
                                         expression=ObservableExpression(expression='st_1[c]', species=[s_1]))
        obs_2 = model.observables.create(id='obs_2',
                                         expression=ObservableExpression(expression='st_2[c] + obs_1', species=[s_2]))
        obs_2.expression.observables = [obs_1]
        rv = model.validate()
        self.assertEqual(rv, None, str(rv))

        model = Model(id='model', version='0.0.1')
        c = model.compartments.create(id='c')
        st_1 = model.species_types.create(id='st_1')
        st_2 = model.species_types.create(id='st_2')
        s_1 = model.species.create(species_type=st_1, compartment=c)
        s_2 = model.species.create(species_type=st_2, compartment=c)
        obs_1 = model.observables.create(id='obs_1',
                                         expression=ObservableExpression(expression='st_1[c] + obs_2', species=[s_1]))
        obs_2 = model.observables.create(id='obs_2',
                                         expression=ObservableExpression(expression='st_2[c] + obs_1', species=[s_2]))
        obs_1.expression.observables = [obs_2]
        obs_2.expression.observables = [obs_1]
        rv = model.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertRegex(str(rv), 'cannot have cyclic depencencies')

        model = Model(id='model', version='0.0.1')
        c = model.compartments.create(id='c')
        st_1 = model.species_types.create(id='st_1')
        st_2 = model.species_types.create(id='st_2')
        s_1 = model.species.create(species_type=st_1, compartment=c)
        s_2 = model.species.create(species_type=st_2, compartment=c)
        obs_1 = model.observables.create(id='obs_1',
                                         expression=ObservableExpression(expression='st_1[c] + obs_1', species=[s_1]))
        obs_2 = model.observables.create(id='obs_2',
                                         expression=ObservableExpression(expression='st_2[c]', species=[s_2]))
        obs_1.expression.observables = [obs_1]
        rv = model.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertRegex(str(rv), 'cannot have cyclic depencencies')

        model = Model(id='model', version='0.0.1')
        c = model.compartments.create(id='c')
        st_1 = model.species_types.create(id='st_1')
        st_2 = model.species_types.create(id='st_2')
        s_1 = model.species.create(species_type=st_1, compartment=c)
        s_2 = model.species.create(species_type=st_2, compartment=c)
        func_1 = model.functions.create(id='func_1',
                                        expression=FunctionExpression(expression='st_1[c] + func_2', species=[s_1]))
        func_2 = model.functions.create(id='func_2',
                                        expression=FunctionExpression(expression='st_2[c] + func_1', species=[s_2]))
        func_1.expression.functions = [func_2]
        func_2.expression.functions = [func_1]
        rv = model.validate()
        self.assertEqual(len(rv.attributes), 1)
        self.assertRegex(str(rv), 'cannot have cyclic depencencies')

    def test_observation_validate(self):
        ev = Observation(id='ev')
        error = ev.validate()
        self.assertEqual(error, None, str(error))

        ev = Observation(id='')
        error = ev.validate()
        self.assertNotEqual(error, None, str(error))

        env = ObservationEnv(temp=1., temp_units=unit_registry.parse_units('celsius'))
        error = env.validate()
        self.assertEqual(error, None, str(error))

        env = ObservationEnv(temp=1.)
        error = env.validate()
        self.assertNotEqual(error, None, str(error))

        env = ObservationEnv(ph=1., ph_units=unit_registry.parse_units('dimensionless'))
        error = env.validate()
        self.assertEqual(error, None, str(error))

        env = ObservationEnv(ph=1.)
        error = env.validate()
        self.assertNotEqual(error, None, str(error))

        env = ObservationEnv(ph=1., ph_units=unit_registry.parse_units('celsius'))
        error = env.validate()
        self.assertNotEqual(error, None, str(error))


class UnitsTestCase(unittest.TestCase):
    def test_species_count(self):
        self.assertEqual(Species.units.choices, Observable.units.choices)

    def test_function_value(self):
        self.assertTrue(hasattr(Function, 'units'))

        objs = {
            SpeciesType: {
                'st_1': SpeciesType(id='st_1'),
                'st_2': SpeciesType(id='st_2'),
            },
            Compartment: {
                'c_1': Compartment(id='c_1'),
                'c_2': Compartment(id='c_2'),
            },
            Species: {
            },
            Observable: {
            },
            Function: {
            },
            Parameter: {
                'p_1': Parameter(id='p_1', value=1.5, units=unit_registry.parse_units('g')),
                'p_2': Parameter(id='p_2', value=2.5, units=unit_registry.parse_units('kg')),
            },
        }
        objs[Species]['st_1[c_1]'] = Species(species_type=objs[SpeciesType]['st_1'],
                                             compartment=objs[Compartment]['c_1'])
        objs[Species]['st_2[c_2]'] = Species(species_type=objs[SpeciesType]['st_2'],
                                             compartment=objs[Compartment]['c_2'])
        objs[Species]['st_1[c_1]'].id = objs[Species]['st_1[c_1]'].gen_id()
        objs[Species]['st_2[c_2]'].id = objs[Species]['st_2[c_2]'].gen_id()
        obs_1 = objs[Observable]['obs_1'] = Observable(id='obs_1')
        obs_1.expression, invalid = ObservableExpression.deserialize('3 * st_1[c_1]', objs)
        obs_2 = objs[Observable]['obs_2'] = Observable(id='obs_2')
        obs_2.expression, _ = ObservableExpression.deserialize('4 * st_2[c_2] + obs_1', objs)

        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('1', objs)[0],
            units=unit_registry.parse_units('dimensionless'))
        rv = function.validate()
        self.assertEqual(rv, None, str(rv))

        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('2', objs)[0],
            units=unit_registry.parse_units('s'))
        rv = function.validate()
        self.assertNotEqual(rv, None, str(rv))
        self.assertRegex(str(rv), 'Units of ".*?" should be ".*?" not "dimensionless"')

        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('2 * st_1[c_1]', objs)[0],
            units=unit_registry.parse_units('molecule'))
        self.assertEqual(function.expression.species, [objs[Species]['st_1[c_1]']])
        self.assertEqual(function.expression.compartments, [])
        rv = function.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(function.expression._parsed_expression.test_eval({Species: 2.}), 4.)
        self.assertEqual(function.expression._parsed_expression.test_eval({Species: 3.}), 6.)

        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('2 * st_1[c_1] / c_1', objs)[0],
            units=unit_registry.parse_units('molecule g^-1'))
        self.assertEqual(function.expression.species, [objs[Species]['st_1[c_1]']])
        self.assertEqual(function.expression.compartments, [objs[Compartment]['c_1']])
        rv = function.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(function.expression._parsed_expression.test_eval(
            {Species: 2., Compartment: 3.}), 4./3.)
        self.assertEqual(function.expression._parsed_expression.test_eval(
            {Species: 3., Compartment: 3.}), 6./3.)

        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('2 * st_1[c_1] / Compartment.c_1', objs)[0],
            units=unit_registry.parse_units('molecule g^-1'))
        self.assertEqual(function.expression.species, [objs[Species]['st_1[c_1]']])
        self.assertEqual(function.expression.compartments, [objs[Compartment]['c_1']])
        rv = function.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(function.expression._parsed_expression.test_eval(
            {Species: 2., Compartment: 3.}), 4./3.)
        self.assertEqual(function.expression._parsed_expression.test_eval(
            {Species: 3., Compartment: 3.}), 6./3.)

        # with parameters
        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('2 * st_1[c_1] / p_1', objs)[0],
            units=unit_registry.parse_units('molecule g^-1'))
        rv = function.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(function.expression._parsed_expression.test_eval({Species: 2., Parameter: {'p_1': 1.5}}), 8./3.)

        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('2 * st_1[c_1] / p_2', objs)[0],
            units=unit_registry.parse_units('molecule kg^-1'))
        rv = function.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(function.expression._parsed_expression.test_eval({Species: 2., Parameter: {'p_2': 2.5}}), 8./5.)

        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('(2 * st_1[c_1] + st_2[c_2]) / p_2', objs)[0],
            units=unit_registry.parse_units('molecule kg^-1'))
        rv = function.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(function.expression._parsed_expression.test_eval({Species: 2., Parameter: {'p_2': 2.5}}), 2.4)

        # inconsistent units
        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('2 * st_1[c_1] + st_2[c_2] / p_2', objs)[0],
            units=unit_registry.parse_units('molecule kg^-1'))
        rv = function.validate()
        self.assertNotEqual(rv, None, str(rv))
        self.assertRegex(str(rv), 'Cannot convert from ')

        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('2 * st_1[c_1] + st_2[c_2] / (p_1 / p_2)', objs)[0],
            units=unit_registry.parse_units('molecule'))
        rv = function.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(function.expression._parsed_expression.test_eval(
            {Species: 2., Parameter: {'p_1': 1.5, 'p_2': 2.5}}, with_units=True),
            (4. + 2. / (1.5/2.5e3)) * unit_registry.parse_units('molecule'))

        function = Function(
            id='func',
            expression=FunctionExpression.deserialize('obs_1', objs)[0],
            units=unit_registry.parse_units('molecule'))
        rv = function.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(function.expression._parsed_expression.test_eval(
            {Species: 2.}), 6.)

        func = Function(
            id='func',
            expression=FunctionExpression.deserialize('obs_2', objs)[0],
            units=unit_registry.parse_units('molecule'))
        rv = func.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(func.expression._parsed_expression.test_eval(
            {Species: 2.}), 14.)

        objs[Function]['func'] = func
        func2 = Function(
            id='func_2',
            units=unit_registry.parse_units('molecule'))
        func2.expression, error = FunctionExpression.deserialize('func + p_1 / p_2 * st_1[c_1]', objs)
        self.assertEqual(error, None, str(error))
        rv = func2.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(func2.expression._parsed_expression.test_eval(
            {Species: 2., Parameter: {'p_1': 1.5, 'p_2': 2.5}}, with_units=True).magnitude, 14. + 1.5/2.5e3 * 2)

        func2 = Function(
            id='func_2',
            units=unit_registry.parse_units('molecule'))
        func2.expression, error = FunctionExpression.deserialize('func() + p_1 / p_2 * st_1[c_1]', objs)
        self.assertEqual(error, None, str(error))
        rv = func2.validate()
        self.assertNotEqual(rv, None, str(rv))

        func2 = Function(
            id='func_2',
            units=unit_registry.parse_units('molecule'))
        func2.expression, error = FunctionExpression.deserialize('func( ) + p_1 / p_2 * st_1[c_1]', objs)
        self.assertEqual(error, None, str(error))
        rv = func2.validate()
        self.assertNotEqual(rv, None, str(rv))

    def test_rate_law_value(self):
        self.assertTrue(hasattr(RateLaw, 'units'))

        objs = {
            SpeciesType: {
                'st_1': SpeciesType(id='st_1'),
            },
            Compartment: {
                'c_1': Compartment(id='c_1'),
            },
            Species: {
            },
            Observable: {
            },
            Function: {
            },
            Parameter: {
                'p_1': Parameter(id='p_1', value=1.5, units=unit_registry.parse_units('molecule^-1 s^-1')),
                'p_2': Parameter(id='p_2', value=1.5, units=unit_registry.parse_units('molecule^-1 g s^-1')),
                'p_3': Parameter(id='p_3', value=1.5, units=unit_registry.parse_units('molecule g^-1 s^-1')),
            },
        }
        objs[Species]['st_1[c_1]'] = Species(species_type=objs[SpeciesType]['st_1'],
                                             compartment=objs[Compartment]['c_1'])
        objs[Species]['st_1[c_1]'].id = objs[Species]['st_1[c_1]'].gen_id()
        obs_1 = objs[Observable]['obs_1'] = Observable(id='obs_1')
        obs_1.expression, invalid = ObservableExpression.deserialize('2 * st_1[c_1]', objs)

        func_1 = objs[Function]['func_1'] = Function(id='func_1', units=unit_registry.parse_units('molecule'))
        func_1.expression, _ = FunctionExpression.deserialize('3 * obs_1', objs)

        rl = RateLaw(id='rxn_1-forward',
                     reaction=Reaction(id='rxn_1'),
                     direction=RateLawDirection.forward,
                     units=unit_registry.parse_units('s^-1'))
        rl.expression, _ = RateLawExpression.deserialize('4 * p_1 * func_1', objs)
        rv = rl.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(rl.expression._parsed_expression.test_eval(
            {Species: 2., Parameter: {'p_1': 1.5}}, with_units=True).magnitude,
            4. * 1.5 * (3. * (2. * 2.)))

        rl = RateLaw(id='rxn_1-forward',
                     reaction=Reaction(id='rxn_1'),
                     direction=RateLawDirection.forward,
                     units=unit_registry.parse_units('s^-1'))
        rl.expression, _ = RateLawExpression.deserialize('4 * func_1', objs)
        rv = rl.validate()
        self.assertNotEqual(rv, None, str(rv))

        l = RateLaw(id='rxn_1-forward',
                    reaction=Reaction(id='rxn_1'),
                    direction=RateLawDirection.forward,
                    units=unit_registry.parse_units('s^-1'))
        rl.expression, _ = RateLawExpression.deserialize('4 * p_2 * st_1[c_1] / c_1', objs)
        self.assertEqual(rl.expression.compartments, [objs[Compartment]['c_1']])
        self.assertEqual(rl.expression.species, [objs[Species]['st_1[c_1]']])
        self.assertEqual(rl.expression.parameters, [objs[Parameter]['p_2']])
        rv = rl.validate()
        self.assertEqual(rv, None, str(rv))

        l = RateLaw(id='rxn_1-forward',
                    reaction=Reaction(id='rxn_1'),
                    direction=RateLawDirection.forward,
                    units=unit_registry.parse_units('s^-1'))
        rl.expression, _ = RateLawExpression.deserialize('4 * p_3 * c_1 / st_1[c_1]', objs)
        self.assertEqual(rl.expression.compartments, [objs[Compartment]['c_1']])
        self.assertEqual(rl.expression.species, [objs[Species]['st_1[c_1]']])
        self.assertEqual(rl.expression.parameters, [objs[Parameter]['p_3']])
        rv = rl.validate()
        self.assertEqual(rv, None, str(rv))

        rl.reaction.rate_units = None
        rv = rl.validate()
        self.assertNotEqual(rv, None, str(rv))

    def test_dfba_obj_value(self):
        self.assertEqual(len(DfbaObjective.units.choices), 1)
        self.assertEqual(len(DfbaObjective.reaction_rate_units.choices), 1)
        self.assertEqual(len(DfbaObjective.coefficient_units.choices), 1)

        self.assertEqual(DfbaObjective.reaction_rate_units.choices[0] * DfbaObjective.coefficient_units.choices[0],
                         DfbaObjective.units.choices[0])

        submodel = Submodel(id='submdl')
        rxn_1 = Reaction(id='rxn_1', submodel=submodel)
        dfba_obj_rxn_1 = DfbaObjReaction(id='dfba_obj_rxn_1', submodel=submodel)
        objs = {
            Reaction: {
                'rxn_1': rxn_1,
            },
            DfbaObjReaction: {
                'dfba_obj_rxn_1': dfba_obj_rxn_1,
            },
        }

        dfba_obj = DfbaObjective(submodel=submodel)
        dfba_obj.id = dfba_obj.gen_id()
        dfba_obj.expression, error = DfbaObjectiveExpression.deserialize('rxn_1 + dfba_obj_rxn_1', objs)
        self.assertEqual(error, None)
        rv = dfba_obj.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertEqual(dfba_obj.expression._parsed_expression.test_eval(
            {Reaction: 2.1, DfbaObjReaction: 3.}, with_units=False),
            (2.1 + 3.))

        dfba_obj.expression._parsed_expression._compiled_expression_with_units = ' + '.join([
            'unit_registry.parse_units("s") * Reaction["rxn_1"]',
            'unit_registry.parse_units("s") * DfbaObjReaction["dfba_obj_rxn_1"]',
        ])
        dfba_obj.expression._parsed_expression._compiled_namespace_with_units['unit_registry'] = unit_registry
        self.assertEqual(dfba_obj.expression._parsed_expression.test_eval(
            {Reaction: 2.1, DfbaObjReaction: 3.}, with_units=True),
            (2.1 + 3.) * unit_registry.parse_units('dimensionless'))

    def test_dfba_obj_specices_value(self):
        time_unit = unit_registry.parse_expression(str(Model.time_units.choices[0]))
        mol_unit = unit_registry.parse_expression('mol')

        coeff_unit = unit_registry.parse_expression('M s^-1')
        cell_size_unit = unit_registry.parse_expression('l')
        self.assertTrue((coeff_unit * cell_size_unit * time_unit).to_base_units() == mol_unit.to_base_units())

        coeff_unit = unit_registry.parse_expression('mol gDCW^-1 s^-1')
        cell_size_unit = unit_registry.parse_expression('gDCW')
        self.assertEqual((coeff_unit * cell_size_unit * time_unit).to_base_units(), mol_unit.to_base_units())

    def test_stop_condition_value(self):
        objs = {
            SpeciesType: {
                'st_1': SpeciesType(id='st_1'),
            },
            Compartment: {
                'c_1': Compartment(id='c_1'),
            },
            Species: {
            },
            Observable: {
            },
            Function: {
            },
            Parameter: {
                'p_1': Parameter(id='p_1', value=1.5, units=unit_registry.parse_units('molecule^-1 s^-1')),
                'p_2': Parameter(id='p_2', value=1.5, units=unit_registry.parse_units('molecule')),
                'p_3': Parameter(id='p_3', value=1.5, units=unit_registry.parse_units('dimensionless')),
                'p_4': Parameter(id='p_4', value=2.5, units=unit_registry.parse_units('dimensionless')),
                'p_5': Parameter(id='p_5', value=1.0, units=unit_registry.parse_units('molecule g^-1')),
            },
        }
        objs[Species]['st_1[c_1]'] = Species(species_type=objs[SpeciesType]['st_1'],
                                             compartment=objs[Compartment]['c_1'])
        objs[Species]['st_1[c_1]'].id = objs[Species]['st_1[c_1]'].gen_id()
        obs_1 = objs[Observable]['obs_1'] = Observable(id='obs_1')
        obs_1.expression, error = ObservableExpression.deserialize('2 * st_1[c_1]', objs)
        self.assertEqual(error, None)
        self.assertEqual(obs_1.expression._parsed_expression.test_eval(values={Species: 1.}), 2. * 1.)

        func_1 = objs[Function]['func_1'] = Function(id='func_1', units=unit_registry.parse_units('molecule'))
        func_1.expression, error = FunctionExpression.deserialize('3 * obs_1', objs)
        self.assertEqual(error, None)
        self.assertEqual(func_1.expression._parsed_expression.test_eval(values={Species: 1.}), 3. * 2. * 1.)

        cond = StopCondition(id='cond_1')
        cond.expression, error = StopConditionExpression.deserialize('func_1 > p_2', objs)
        self.assertEqual(error, None)
        rv = cond.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertTrue(cond.expression._parsed_expression.test_eval(values={Species: 1.}, with_units=True))

        cond = StopCondition(id='cond_2')
        cond.expression, error = StopConditionExpression.deserialize('func_1 > p_3', objs)
        self.assertEqual(error, None)
        rv = cond.validate()
        self.assertNotEqual(rv, None, str(rv))

        cond = StopCondition(id='cond_3')
        cond.expression, error = StopConditionExpression.deserialize('p_3 > p_4', objs)
        self.assertEqual(error, None)
        rv = cond.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertFalse(cond.expression._parsed_expression.test_eval(
            values={Parameter: {'p_3': 1.5, 'p_4': 2.5}}, with_units=True))

        cond = StopCondition(id='cond_4')
        cond.expression, error = StopConditionExpression.deserialize('st_1[c_1] / c_1 > p_5', objs)
        self.assertEqual(error, None, str(error))
        rv = cond.validate()
        self.assertEqual(rv, None, str(rv))
        self.assertTrue(cond.expression._parsed_expression.test_eval(
            values={Species: 2., Compartment: 1., Parameter: {'p_5': 1.}}, with_units=True))
        self.assertFalse(cond.expression._parsed_expression.test_eval(
            values={Species: 0.5, Compartment: 1., Parameter: {'p_5': 1.}}, with_units=True))

    def test_parameter_value(self):
        self.assertTrue(hasattr(Parameter, 'units'))


class ValidatorTestCase(unittest.TestCase):
    def test(self):
        model = Model(id='model', name='test model', version='0.0.1', wc_lang_version='0.0.1')
        self.assertEqual(Validator().run(model), None)

        model.id = ''
        self.assertNotEqual(Validator().run(model), None)
