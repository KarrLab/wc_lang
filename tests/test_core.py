""" Tests of core

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

import os
import pytest
import unittest
import warnings
import wc_lang
from wc_lang.core import (Model, Taxon, TaxonRank, Submodel, ObjectiveFunction,
                          Reaction, SpeciesType, SpeciesTypeType, Species, Observable, Compartment,
                          SpeciesCoefficient, ObservableCoefficient, Parameter, Reference, ReferenceType,
                          DatabaseReference,
                          RateLaw, RateLawEquation, SubmodelAlgorithm, Concentration, BiomassComponent,
                          BiomassReaction, StopCondition,
                          OneToOneSpeciesAttribute, ReactionParticipantAttribute, RateLawEquationAttribute,
                          InvalidObject, EXTRACELLULAR_COMPARTMENT_ID)
from wc_lang.prepare import PrepareModel
from wc_lang.io import Reader
from libsbml import (SBMLNamespaces, SBMLDocument, readSBMLFromString)
import libsbml
from wc_lang.sbml.util import (wrap_libsbml, LibSBMLError, init_sbml_model,
                               create_sbml_doc_w_fbc, SBML_LEVEL, SBML_VERSION, get_SBML_compatibility_method)


class TestCore(unittest.TestCase):

    def setUp(self):
        Reaction.objects.reset()
        BiomassReaction.objects.reset()

        self.model = mdl = Model(id='model', name='test model', version='0.0.1', wc_lang_version='0.0.1')

        mdl.taxon = Taxon(id='taxon', name='test taxon', rank=TaxonRank.species)

        self.comp_0 = comp_0 = mdl.compartments.create(id='comp_0', name='compartment 0',
                                                       initial_volume=1.25)
        self.comp_1 = comp_1 = mdl.compartments.create(id='comp_1', name='compartment 1',
                                                       initial_volume=2.5)
        self.compartments = compartments = [comp_0, comp_1]

        self.species_types = species_types = []
        self.species = species = []
        self.concentrations = concentrations = []
        for i in range(8):
            spec_type = mdl.species_types.create(
                id='spec_type_{}'.format(i),
                name='species type {}'.format(i),
                type=SpeciesTypeType.metabolite,
                structure='C' * i + 'H' * (i + 1),
                empirical_formula='C{}H{}'.format(i, i + 1),
                molecular_weight=12 * (i + 1),
                charge=i + 1)
            species_types.append(spec_type)

            if i != 3:
                spec = Species(species_type=spec_type, compartment=comp_0)
            else:
                spec = Species(species_type=spec_type, compartment=comp_1)
            species.append(spec)

            conc = Concentration(species=spec, value=3 * i)
            concentrations.append(conc)

        self.biomass_reaction = biomass_reaction = BiomassReaction(
            id='biomass_reaction_1',
            name='biomass reaction',
            compartment=comp_0,
            comments="Nobody will ever deprive the American people of the right to vote except the "
            "American people themselves")
        BiomassReaction.get_manager().insert_all_new()

        biomass_components = []
        for i in range(2):
            biomass_components.append(
                biomass_reaction.biomass_components.create(
                    id='biomass_comp_{}'.format(i + 1),
                    coefficient=2 * (float(i) - 0.5),     # create a reactant and a product
                    species_type=species_types[i]))
        self.biomass_components = biomass_components

        self.submdl_0 = submdl_0 = mdl.submodels.create(
            id='submodel_0', name='submodel 0', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_1 = submdl_1 = mdl.submodels.create(
            id='submodel_1', name='submodel 1', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_2 = submdl_2 = mdl.submodels.create(
            id='submodel_2', name='submodel 2', algorithm=SubmodelAlgorithm.dfba, compartment=comp_0,
            biomass_reaction=biomass_reaction)
        self.submodels = submodels = [submdl_0, submdl_1, submdl_2]

        self.rxn_0 = rxn_0 = submdl_0.reactions.create(id='rxn_0', name='reaction 0')
        rxn_0.participants.create(species=species[0], coefficient=-2)
        rxn_0.participants.create(species=species[1], coefficient=-3.5)
        rxn_0.participants.create(species=species[2], coefficient=1)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[5].serialize()),
            modifiers=species[5:6])
        rate_law_0 = rxn_0.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_1 = rxn_1 = submdl_1.reactions.create(id='rxn_1', name='reaction 1')
        rxn_1.participants.create(species=species[0], coefficient=-2)
        rxn_1.participants.create(species=species[1], coefficient=-3)
        rxn_1.participants.create(species=species[3], coefficient=2)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[6].serialize()),
            modifiers=species[6:7])
        rate_law_1 = rxn_1.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_2 = rxn_2 = submdl_2.reactions.create(id='rxn_2', name='reaction 2')
        rxn_2.participants.create(species=species[0], coefficient=-2)
        rxn_2.participants.create(species=species[1], coefficient=-3)
        rxn_2.participants.create(species=species[4], coefficient=1)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[7].serialize()),
            modifiers=species[7:8])
        rate_law_2 = rxn_2.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        Reaction.get_manager().insert_all_new()

        self.reactions = [rxn_0, rxn_1, rxn_2]
        self.rate_laws = [rate_law_0, rate_law_1, rate_law_2]

        self.objective_function = of = biomass_reaction.objective_functions.create()
        of.submodels.append(submdl_2)
        of.reactions.append(rxn_1)
        of.reactions.append(rxn_2)

        self.parameters = parameters = []
        self.references = references = []
        self.database_references = database_references = []
        for i in range(3):
            param = mdl.parameters.create(
                id='param_{}'.format(i), name='parameter {}'.format(i),
                value=i * 4, units='dimensionless')
            param.submodels = submodels[i:i + 1]
            parameters.append(param)

            ref = param.references.create(
                id='ref_{}'.format(i), name='reference {}'.format(i),
                type=ReferenceType.misc)
            references.append(ref)

            x_ref = ref.database_references.create(database='x', id='y' * (i + 1),
                                                   url='http://x.com/{}'.format('y' * (i + 1)))
            database_references.append(x_ref)

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
            self.assertEqual(submodel.database_references, [])
            self.assertEqual(submodel.references, [])

        # compartment
        self.assertEqual(set(self.compartments[0].species), set(self.species[0:3] + self.species[4:]))
        self.assertEqual(self.compartments[1].species, self.species[3:4])

        for compartment in self.compartments:
            self.assertEqual(compartment.database_references, [])
            self.assertEqual(compartment.references, [])

        # species type
        for species_type, species in zip(self.species_types, self.species):
            self.assertEqual(species_type.species, [species])

        for species_type in self.species_types:
            self.assertEqual(species_type.database_references, [])
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

        self.assertEqual(len(self.species[0].rate_law_equations), 0)
        self.assertEqual(len(self.species[1].rate_law_equations), 0)
        self.assertEqual(len(self.species[2].rate_law_equations), 0)
        self.assertEqual(len(self.species[3].rate_law_equations), 0)
        self.assertEqual(len(self.species[4].rate_law_equations), 0)
        self.assertEqual(len(self.species[5].rate_law_equations), 1)
        self.assertEqual(len(self.species[6].rate_law_equations), 1)
        self.assertEqual(len(self.species[7].rate_law_equations), 1)

        # reaction
        for reaction, submodel in zip(self.reactions, self.submodels):
            self.assertEqual(reaction.submodel, submodel)

        self.assertEqual(set(x.species for x in self.reactions[0].participants),
                         set([self.species[0], self.species[1], self.species[2]]))
        self.assertEqual(set(x.species for x in self.reactions[1].participants),
                         set([self.species[0], self.species[1], self.species[3]]))
        self.assertEqual(set(x.species for x in self.reactions[2].participants),
                         set([self.species[0], self.species[1], self.species[4]]))

        self.assertEqual(self.reactions[0].rate_laws[0].equation.modifiers, self.species[5:6])
        self.assertEqual(self.reactions[1].rate_laws[0].equation.modifiers, self.species[6:7])
        self.assertEqual(self.reactions[2].rate_laws[0].equation.modifiers, self.species[7:8])

        for reaction in self.reactions:
            self.assertEqual(reaction.database_references, [])
            self.assertEqual(reaction.references, [])
            self.assertEqual(len(reaction.rate_laws), 1)

        # biomass components
        for i in range(len(self.biomass_components)):
            # submodels
            self.assertEqual(self.biomass_components[i].biomass_reaction, self.biomass_reaction)
            # self.assertEqual(self.biomass_reaction.submodels[0], self.submodels[2])
            # species types
            self.assertEqual(self.biomass_components[i].species_type, self.species_types[i])
            self.assertEqual(self.biomass_components[i], self.species_types[i].biomass_components[0])

        # parameters
        for reference, parameter in zip(self.references, self.parameters):
            self.assertEqual(parameter.references, [reference])

        for parameter in self.parameters:
            self.assertEqual(parameter.model, mdl)

        # references
        for reference, parameter in zip(self.references, self.parameters):
            self.assertEqual(reference.parameters, [parameter])
            self.assertEqual(parameter.references, [reference])

        for reference, database_reference in zip(self.references, self.database_references):
            self.assertEqual(reference.database_references, [database_reference])
            self.assertEqual(database_reference.reference, reference)

        # reaction participant
        for species in self.species[0:5]:
            self.assertEqual(set(x.species for x in species.species_coefficients), set([species]))

        for reaction in self.reactions:
            for part in reaction.participants:
                self.assertIn(reaction, part.reactions)
            self.assertEqual(set(x.reaction for x in reaction.rate_laws), set([reaction]))

        # database references
        for reference, database_reference in zip(self.references, self.database_references):
            self.assertEqual(reference.database_references, [database_reference])
            self.assertEqual(database_reference.reference, reference)

    def test_taxon_rank_class(self):
        self.assertEqual(TaxonRank['class'], TaxonRank['classis'])
        self.assertEqual(TaxonRank.__getattr__('class'), TaxonRank['classis'])

    def test_model_get_species(self):
        self.assertEqual(set(self.model.get_species()), set(self.species))

    def test_submodel_get_species(self):
        species = self.species
        self.assertEqual(set(self.submdl_0.get_species()), set([
            species[0], species[1], species[2], species[5],
        ]))
        self.assertEqual(set(self.submdl_1.get_species()), set([
            species[0], species[1], species[3], species[6],
        ]))
        self.assertEqual(set(self.submdl_2.get_species()), set([
            species[0], species[1], species[4], species[7],
        ]))

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

        self.assertEqual(set(mdl.get_concentrations()), set(self.concentrations))
        self.assertEqual(set(mdl.get_concentrations(__type=Concentration)), set(self.concentrations))
        self.assertEqual(set(mdl.get_concentrations(__type=Submodel)), set())

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

        self.assertNotEqual(set(mdl.get_biomass_reactions(__type=BiomassReaction)), set())
        self.assertEqual(set(mdl.get_biomass_reactions(__type=Reaction)), set())

        self.assertEqual(set(self.objective_function.get_products()), set([
            self.species[3],
            self.species[4],
            Species.get([Species.gen_id(self.species_types[1], self.biomass_reaction.compartment)], mdl.get_species())[0],
        ]))
        self.assertEqual(set(self.objective_function.get_products(__type=Species)), set([
            self.species[3],
            self.species[4],
            Species.get([Species.gen_id(self.species_types[1], self.biomass_reaction.compartment)], mdl.get_species())[0],
        ]))
        self.assertEqual(set(self.objective_function.get_products(__type=Reaction)), set())

    def test_get_component(self):
        model = self.model

        self.assertEqual(model.get_component('compartment', 'comp_0'), self.comp_0)
        self.assertEqual(model.get_component('species_type', 'spec_type_1'), self.species_types[1])
        self.assertEqual(model.get_component('submodel', 'submodel_1'), self.submdl_1)
        self.assertEqual(model.get_component('reaction', 'rxn_1'), self.rxn_1)
        self.assertEqual(model.get_component('parameter', 'param_2'), self.parameters[2])
        self.assertEqual(model.get_component('reference', 'ref_1'), self.references[1])
        self.assertEqual(model.get_component('reaction', 'rxn_3'), None)

        with self.assertRaisesRegexp(ValueError, ' not one of '):
            model.get_component('undefined', 'rxn_3')

    def test_species_type_is_carbon_containing(self):
        self.assertFalse(self.species_types[0].is_carbon_containing())
        self.assertTrue(self.species_types[1].is_carbon_containing())

    def test_species_serialize(self):
        self.assertEqual(self.species[0].serialize(), 'spec_type_0[comp_0]')
        self.assertEqual(self.species[1].serialize(), 'spec_type_1[comp_0]')
        self.assertEqual(self.species[2].serialize(), 'spec_type_2[comp_0]')
        self.assertEqual(self.species[3].serialize(), 'spec_type_3[comp_1]')

    def test_species_gen_id(self):
        self.assertEqual(Species.gen_id(self.species[3].species_type, self.species[3].compartment),
                         'spec_type_3[comp_1]')
        self.assertEqual(
            Species.gen_id(self.species[3].species_type.id, self.species[3].compartment.id),
            'spec_type_3[comp_1]')
        with self.assertRaises(ValueError) as context:
            Species.gen_id(self.species[3].species_type.id, self.species[3].compartment)
        self.assertIn('gen_id: incorrect parameter types', str(context.exception))

    def test_species_get(self):
        self.assertEqual(Species.get([], self.species), [])
        self.assertEqual(Species.get(['X'], self.species), [None])
        self.assertEqual(Species.get(['spec_type_0[comp_0]'], self.species), [self.species[0]])
        ids = ["spec_type_{}[comp_0]".format(i) for i in range(4, 8)]
        self.assertEqual(Species.get(ids, self.species), self.species[4:])
        ids.append('X')
        self.assertEqual(Species.get(ids, self.species), self.species[4:] + [None])

    def test_species_deserialize(self):
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
        }

        attr = Concentration.Meta.attributes['species']

        val = 'spec_0[c_0]'
        species0, error = Species.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(species0.serialize(), val)
        self.assertEqual(set(objs[Species].values()), set([species0]))

        val = 'spec_2[c_1]'
        species1, error = Species.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(species1.serialize(), val)
        self.assertEqual(set(objs[Species].values()), set([species0, species1]))

        val = 'spec_2[c_3]'
        species2, error = Species.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(species2, None)
        self.assertEqual(set(objs[Species].values()), set([species0, species1]))

        val = 'spec_2'
        species3, error = Species.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(species3, None)
        self.assertEqual(set(objs[Species].values()), set([species0, species1]))

        val = '[c_3]'
        species4, error = Species.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(species4, None)
        self.assertEqual(set(objs[Species].values()), set([species0, species1]))

    def test_observable_species_serialize(self):
        st_a = SpeciesType(id='a')
        st_b = SpeciesType(id='bb')
        st_c = SpeciesType(id='ccc')
        st_d = SpeciesType(id='dddd')
        c_a = Compartment(id='a')
        c_b = Compartment(id='bb')
        c_c = Compartment(id='ccc')
        c_d = Compartment(id='dddd')
        s_a = Species(species_type=st_a, compartment=c_a)
        s_b = Species(species_type=st_b, compartment=c_b)
        s_c = Species(species_type=st_c, compartment=c_c)
        s_d = Species(species_type=st_d, compartment=c_d)
        sc_a = SpeciesCoefficient(species=s_a, coefficient=2.)
        sc_b = SpeciesCoefficient(species=s_b, coefficient=3.)
        sc_c = SpeciesCoefficient(species=s_c, coefficient=4.)
        sc_d = SpeciesCoefficient(species=s_d, coefficient=1.)
        obs = Observable()
        obs.species.append(sc_a)
        obs.species.append(sc_b)
        obs.species.append(sc_c)
        obs.species.append(sc_d)

        objs = {
            SpeciesType: {
                st_a.serialize(): st_a,
                st_b.serialize(): st_b,
                st_c.serialize(): st_c,
                st_d.serialize(): st_d,
            },
            Compartment: {
                c_a.serialize(): c_a,
                c_b.serialize(): c_b,
                c_c.serialize(): c_c,
                c_d.serialize(): c_d,
            },
            Species: {
                s_a.serialize(): s_a,
                s_b.serialize(): s_b,
            },
            SpeciesCoefficient: {
                sc_a.serialize(): sc_a,
            }
        }

        attr = Observable.Meta.attributes['species']

        self.assertEqual(sc_a.serialize(), '(2) a[a]')
        self.assertEqual(sc_b.serialize(), '(3) bb[bb]')
        self.assertEqual(sc_d.serialize(), 'dddd[dddd]')
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(2) a[a]', objs)[0].species.species_type, st_a)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(2) a[a]', objs)[0].species.compartment, c_a)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(2) a[a]', objs)[0].coefficient, 2.)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(2) a[a]', objs)[0].species, s_a)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(2) a[a]', objs)[0], sc_a)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(3) bb[bb]', objs)[0].species.species_type, st_b)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(3) bb[bb]', objs)[0].species.compartment, c_b)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(3) bb[bb]', objs)[0].coefficient, 3.)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(3) bb[bb]', objs)[0].species, s_b)
        self.assertNotEqual(SpeciesCoefficient.deserialize(attr, '(3) bb[bb]', objs)[0], sc_b)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(4) ccc[ccc]', objs)[0].species.species_type, st_c)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(4) ccc[ccc]', objs)[0].species.compartment, c_c)
        self.assertEqual(SpeciesCoefficient.deserialize(attr, '(4) ccc[ccc]', objs)[0].coefficient, 4.)
        self.assertNotEqual(SpeciesCoefficient.deserialize(attr, '(4) ccc[ccc]', objs)[0].species, s_c)
        self.assertNotEqual(SpeciesCoefficient.deserialize(attr, '(4) ccc[ccc]', objs)[0], sc_c)

        self.assertEqual(attr.serialize(obs.species), '(2) a[a] + (3) bb[bb] + (4) ccc[ccc] + dddd[dddd]')
        result = attr.deserialize('(2) a[a] + (3) bb[bb] + (4) ccc[ccc] + dddd[dddd]', objs)[0]

        self.assertEqual(result[0], sc_a)
        self.assertEqual(result[1].species, s_b)
        self.assertEqual(result[1].coefficient, 3.)
        self.assertEqual(result[2].species.species_type, st_c)
        self.assertEqual(result[2].species.compartment, c_c)
        self.assertEqual(result[2].coefficient, 4.)
        self.assertEqual(result[3].species.species_type, st_d)
        self.assertEqual(result[3].species.compartment, c_d)
        self.assertEqual(result[3].coefficient, 1.)

        objs = {
            SpeciesType: {
                st_a.id: st_a,
                st_b.id: st_b,
                st_c.id: st_c,
            },
            Compartment: {
                c_a.id: c_a,
                c_b.id: c_b,
                c_c.id: c_c,
            },
        }
        self.assertEqual(result[0].species.species_type, st_a)
        self.assertEqual(result[0].species.compartment, c_a)
        self.assertEqual(result[0].coefficient, 2.)
        self.assertEqual(result[1].species.species_type, st_b)
        self.assertEqual(result[1].species.compartment, c_b)
        self.assertEqual(result[1].coefficient, 3.)
        self.assertEqual(result[2].species.species_type, st_c)
        self.assertEqual(result[2].species.compartment, c_c)
        self.assertEqual(result[2].coefficient, 4.)

        self.assertEqual(attr.serialize([]), '')

        # test deserialize error handling
        self.assertNotEqual(result, None)
        self.assertEqual(attr.deserialize('(2) a[a] - (3) bb[bb] + (4) ccc[ccc]', objs)[0], None)
        self.assertEqual(attr.deserialize('(2) aa[a] + (3) bb[bb] + (4) ccc[ccc]', objs)[0], None)
        self.assertEqual(attr.deserialize('(2) a[aa] + (3) bb[bb] + (4) ccc[ccc]', objs)[0], None)

    def test_observable_observable_serialize(self):
        obs_a = Observable(id='obs_a')
        obs_b = Observable(id='obs_b')
        obs_c = Observable(id='obs_c')

        obs_coeff_a = ObservableCoefficient(observable=obs_a, coefficient=2.)
        obs_coeff_b = ObservableCoefficient(observable=obs_b, coefficient=3.5)
        obs_coeff_c = ObservableCoefficient(observable=obs_c, coefficient=1.)

        self.assertEqual(obs_coeff_a.serialize(), '(2) obs_a')
        self.assertEqual(obs_coeff_b.serialize(), '(3.500000e+00) obs_b')
        self.assertEqual(obs_coeff_c.serialize(), 'obs_c')

        obs_ab = Observable(id='obs_ab')
        obs_ab.observables.append(obs_coeff_a)
        obs_ab.observables.append(obs_coeff_b)

        objs = {
            Observable: {
                obs_a.id: obs_a,
                obs_b.id: obs_b,
            },
            ObservableCoefficient: {
                obs_coeff_a.serialize(): obs_coeff_a,
            },
        }

        attr = Observable.Meta.attributes['observables']

        self.assertEqual(ObservableCoefficient.deserialize(attr, '(2) obs_a', objs)[0], obs_coeff_a)
        self.assertEqual(ObservableCoefficient.deserialize(attr, '(3.5) obs_b', objs)[0].observable, obs_b)
        self.assertEqual(ObservableCoefficient.deserialize(attr, '(3.5) obs_b', objs)[0].coefficient, 3.5)

        objs = {
            Observable: {
                obs_a.id: obs_a,
                obs_b.id: obs_b,
            },
        }
        self.assertEqual(ObservableCoefficient.deserialize(attr, '(3.5) obs_d', objs)[0], None)
        self.assertEqual(ObservableCoefficient.deserialize(attr, '(3.5) obs_b', objs)[0].observable, obs_b)
        self.assertEqual(ObservableCoefficient.deserialize(attr, '(3.5) obs_b', objs)[0].coefficient, 3.5)
        self.assertEqual(ObservableCoefficient.deserialize(attr, '(3.5) obs_b[a]', objs)[0], None)

        objs = {
            Observable: {
                obs_a.id: obs_a,
                obs_b.id: obs_b,
            },
            ObservableCoefficient: {
                obs_coeff_a.serialize(): obs_coeff_a,
            },
        }
        self.assertEqual(attr.serialize(obs_ab.observables), '(2) obs_a + (3.500000e+00) obs_b')
        self.assertEqual(attr.deserialize('(2) obs_a + (3.5) obs_b', objs)[0][0], obs_coeff_a)
        self.assertEqual(attr.deserialize('(2) obs_a + (3.5) obs_b', objs)[0][1].observable, obs_b)
        self.assertEqual(attr.deserialize('(2) obs_a + (3.5) obs_b', objs)[0][1].coefficient, 3.5)

        self.assertEqual(attr.deserialize('(2) obs_a - (3.5) obs_b', objs)[0], None)
        self.assertEqual(attr.deserialize('(2) obs_d + (3.5) obs_b', objs)[0], None)

    def test_concentration_serialize(self):
        self.assertEqual(self.concentrations[0].serialize(), 'spec_type_0[comp_0]')
        self.assertEqual(self.concentrations[1].serialize(), 'spec_type_1[comp_0]')
        self.assertEqual(self.concentrations[2].serialize(), 'spec_type_2[comp_0]')
        self.assertEqual(self.concentrations[3].serialize(), 'spec_type_3[comp_1]')

    def test_reaction_participant_serialize(self):
        self.assertEqual(set([part.serialize() for part in self.rxn_0.participants]), set([
            '(-2) spec_type_0[comp_0]', '(-3.500000e+00) spec_type_1[comp_0]', 'spec_type_2[comp_0]'
        ]))

    def test_reaction_participant_deserialize(self):
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
        }

        attr = Reaction.Meta.attributes['participants']

        val = 'spec_0[c_0]'
        part0, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part0.coefficient, 1)
        self.assertEqual(part0.species.serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[SpeciesCoefficient].values()), set([part0]))
        self.assertEqual(set(objs[Species].values()), set([part0.species]))

        val = '(2) spec_0[c_0]'
        part1, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part1.coefficient, 2)
        self.assertEqual(part1.species.serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[SpeciesCoefficient].values()), set([part0, part1]))
        self.assertEqual(set(objs[Species].values()), set([part0.species, part1.species]))

        val = '(2.) spec_0[c_1]'
        part2, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part2.coefficient, 2)
        self.assertEqual(part2.species.serialize(), 'spec_0[c_1]')
        self.assertEqual(set(objs[SpeciesCoefficient].values()), set([part0, part1, part2]))
        self.assertEqual(set(objs[Species].values()), set([part0.species, part1.species, part2.species]))

        val = '(2.5) spec_0[c_0]'
        part3, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part3.coefficient, 2.5)
        self.assertEqual(part3.species.serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[SpeciesCoefficient].values()), set([part0, part1, part2, part3]))
        self.assertEqual(set(objs[Species].values()), set([part0.species, part1.species, part2.species, part3.species]))

        val = '(.5) spec_0[c_0]'
        part4, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part4.coefficient, 0.5)
        self.assertEqual(part4.species.serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[SpeciesCoefficient].values()), set([part0, part1, part2, part3, part4]))
        self.assertEqual(set(objs[Species].values()), set(
            [part0.species, part1.species, part2.species, part3.species, part4.species]))

        val = '(1) spec_1'
        part5, error = SpeciesCoefficient.deserialize(attr, val, objs, compartment=objs[Compartment]['c_0'])
        self.assertEqual(error, None)
        self.assertEqual(part5.coefficient, 1)
        self.assertEqual(part5.species.serialize(), 'spec_1[c_0]')
        self.assertEqual(set(objs[SpeciesCoefficient].values()), set([part0, part1, part2, part3, part4, part5]))
        self.assertEqual(set(objs[Species].values()), set(
            [part0.species, part1.species, part2.species, part3.species, part4.species, part5.species]))

        # negative examples
        val = '(-1) spec_0[c_0]'
        part6, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        val = '(1) spec_0'
        part6, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        val = '(1.1.) spec_0[c_0]'
        part6, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        val = ' spec_0[c_0]'
        part6, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        val = ' spec_3[c_0]'
        part6, error = SpeciesCoefficient.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        self.assertEqual(set(objs[SpeciesCoefficient].values()), set([part0, part1, part2, part3, part4, part5]))
        self.assertEqual(set(objs[Species].values()), set(
            [part0.species, part1.species, part2.species, part3.species, part4.species, part5.species]))

        val = '(1) spec_3'
        part, error = SpeciesCoefficient.deserialize(attr, val, objs, compartment=objs[Compartment]['c_0'])
        self.assertNotEqual(error, None)
        self.assertEqual(part, None)

        val = '(2) spec_0'
        objs[SpeciesCoefficient] = {
            '(2) spec_0[c_0]': SpeciesCoefficient(
                species=Species(species_type=objs[SpeciesType]['spec_0'], compartment=objs[Compartment]['c_0']),
                coefficient=2)
        }
        part, error = SpeciesCoefficient.deserialize(attr, val, objs, compartment=objs[Compartment]['c_0'])
        self.assertEqual(error, None)
        self.assertEqual(part, objs[SpeciesCoefficient]['(2) spec_0[c_0]'])

    def test_rate_law_serialize(self):
        self.assertEqual(self.rate_laws[0].serialize(), 'rxn_0.forward')
        self.assertEqual(self.rate_laws[1].serialize(), 'rxn_1.forward')
        self.assertEqual(self.rate_laws[2].serialize(), 'rxn_2.forward')

    def test_rate_law_equation_serialize(self):
        self.assertEqual(self.rate_laws[0].equation.serialize(),
                         'k_cat * {0} / (k_m + {0})'.format(self.species[5].serialize()))
        self.assertEqual(self.rate_laws[1].equation.serialize(),
                         'k_cat * {0} / (k_m + {0})'.format(self.species[6].serialize()))
        self.assertEqual(self.rate_laws[2].equation.serialize(),
                         'k_cat * {0} / (k_m + {0})'.format(self.species[7].serialize()))

    def test_rate_law_equation_deserialize(self):
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
        }

        expression = 'k_cat * spec_0[c_0]'
        attr = RateLaw.Meta.attributes['equation']
        equation1, error = RateLawEquation.deserialize(attr, expression, objs)
        self.assertEqual(error, None)
        self.assertEqual(equation1.expression, expression)
        self.assertEqual(equation1.modifiers[0].serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[RateLawEquation].values()), set([equation1]))

        expression = 'k_cat * spec_0[c_1] / (k_m + spec_2[c_1])'
        attr = RateLaw.Meta.attributes['equation']
        equation2, error = RateLawEquation.deserialize(attr, expression, objs)
        self.assertEqual(error, None)
        self.assertEqual(equation2.expression, expression)
        self.assertEqual(set([x.serialize() for x in equation2.modifiers]), set(['spec_0[c_1]', 'spec_2[c_1]']))
        self.assertEqual(set(objs[RateLawEquation].values()), set([equation1, equation2]))

        expression = 'k_cat * spec_0[c_3] / (k_m + spec_1[c_1])'
        attr = RateLaw.Meta.attributes['equation']
        equation, error = RateLawEquation.deserialize(attr, expression, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(equation, None)
        self.assertEqual(set(objs[RateLawEquation].values()), set([equation1, equation2]))

        expression = 'k_cat * spec_3[c_0] / (k_m + spec_1[c_1])'
        attr = RateLaw.Meta.attributes['equation']
        equation, error = RateLawEquation.deserialize(attr, expression, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(equation, None)
        self.assertEqual(set(objs[RateLawEquation].values()), set([equation1, equation2]))

        # exception
        attr = RateLaw.Meta.attributes['equation']
        equation, error = RateLawEquation.deserialize(attr, 2, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(equation, None)

    def test_rate_law_validate(self):
        species_types = [
            SpeciesType(id='spec_0'),
            SpeciesType(id='spec_1'),
        ]
        compartments = [
            Compartment(id='c_0'),
            Compartment(id='c_1'),
        ]

        # unknown specie error
        expression = 'spec_x[c_0]'
        equation = RateLawEquation(
            expression=expression,
            modifiers=[
                Species(species_type=species_types[0], compartment=compartments[0])
            ])
        rate_law = RateLaw(
            equation=equation,
        )
        self.assertNotEqual(rate_law.validate(), None)

        # Name error
        expression = 'not_k_cat * spec_0[c_0]'
        equation = RateLawEquation(
            expression=expression,
            modifiers=[
                Species(species_type=species_types[0], compartment=compartments[0])
            ])
        rate_law = RateLaw(
            equation=equation
        )
        self.assertNotEqual(rate_law.validate(), None)

        # syntax error
        expression = '* spec_0[c_0]'
        equation = RateLawEquation(
            expression=expression,
            modifiers=[
                Species(species_type=species_types[0], compartment=compartments[0])
            ])
        rate_law = RateLaw(
            equation=equation
        )
        self.assertNotEqual(rate_law.validate(), None)

        # No error
        expression = 'k_cat * spec_0[c_0]'
        equation = RateLawEquation(
            expression=expression,
            modifiers=[
                Species(species_type=species_types[0], compartment=compartments[0])
            ])
        rate_law = RateLaw(
            k_cat=2,
            k_m=1,
            equation=equation
        )
        self.assertEqual(rate_law.validate(), None)

    def test_rate_law_equation_validate(self):
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

        expression = 'spec_0[c_0]'
        equation = RateLawEquation(
            expression=expression,
            modifiers=[
                Species(species_type=species_types[0], compartment=compartments[0])
            ])
        self.assertEqual(equation.validate(), None)

        expression = 'spec_0[c_0] * spec_1[c_2]'
        equation = RateLawEquation(
            expression=expression,
            modifiers=[
                Species(species_type=species_types[0], compartment=compartments[0]),
                Species(species_type=species_types[1], compartment=compartments[2]),
            ])
        self.assertEqual(equation.validate(), None)

        expression = 'spec_0[c_0] * spec_1[c_2]'
        equation = RateLawEquation(
            expression=expression,
            modifiers=[
                Species(species_type=species_types[0], compartment=compartments[0]),
                Species(species_type=species_types[1], compartment=compartments[1]),
                Species(species_type=species_types[1], compartment=compartments[2]),
            ])
        self.assertNotEqual(equation.validate(), None)

        expression = 'spec_0[c_0] * spec_1[c_2]'
        equation = RateLawEquation(
            expression=expression,
            modifiers=[
                Species(species_type=species_types[0], compartment=compartments[0]),
            ])
        self.assertNotEqual(equation.validate(), None)

    def test_rate_law_modifiers(self):
        self.assertEqual(self.rxn_0.rate_laws[0].equation.modifiers, self.species[5:6])
        self.assertEqual(self.rxn_1.rate_laws[0].equation.modifiers, self.species[6:7])
        self.assertEqual(self.rxn_2.rate_laws[0].equation.modifiers, self.species[7:8])

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
            Parameter(id='a', submodels=[submodel]),
            Parameter(id='b', submodels=[submodel]),
        ]
        self.assertEqual(Parameter.validate_unique(params), None)

        submodel = Submodel()
        params = [
            Parameter(id='a', submodels=[submodel]),
            Parameter(id='a', submodels=[submodel]),
        ]
        self.assertNotEqual(Parameter.validate_unique(params), None)

        model = Model()
        submodel = Submodel()
        params = [
            Parameter(id='a', model=model),
            Parameter(id='a', submodels=[submodel]),
        ]
        self.assertEqual(Parameter.validate_unique(params), None)

        params = [
            Parameter(id='a', submodels=[Submodel(id='a')]),
            Parameter(id='a', submodels=[Submodel(id='b')]),
        ]
        self.assertEqual(Parameter.validate_unique(params), None)

        params = [
            Parameter(id='a', submodels=[Submodel()]),
            Parameter(id='a', submodels=[Submodel()]),
        ]
        self.assertNotEqual(Parameter.validate_unique(params), None)

    def test_database_reference_serialize(self):
        self.assertEqual(self.database_references[0].serialize(), '{}: {}'.format('x', 'y'))
        self.assertEqual(self.database_references[1].serialize(), '{}: {}'.format('x', 'yy'))
        self.assertEqual(self.database_references[2].serialize(), '{}: {}'.format('x', 'yyy'))

    def test_OneToOneSpeciesAttribute_serialize(self):
        attr = OneToOneSpeciesAttribute()
        self.assertEqual(attr.serialize(self.species[0]), 'spec_type_0[comp_0]')

    def test_OneToOneSpeciesAttribute_deserialize(self):
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
        }

        val = 'spec_0[c_0]'
        attr = OneToOneSpeciesAttribute()
        species0, error = attr.deserialize(val, objs)
        self.assertEqual(error, None)
        self.assertEqual(species0.serialize(), val)
        self.assertEqual(list(objs[Species].values()), [species0])

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
        }

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
        self.assertEqual(len(objs[Species]), 3)
        self.assertEqual(set(objs[Species].values()), set([p.species for p in parts1]))

        parts2, error = attr.deserialize(
            '(2) spec_0[c_0] + (3) spec_1[c_0] ==> (2) spec_2[c_1]', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts2]), set(
            ['(-2) spec_0[c_0]', '(-3) spec_1[c_0]', '(2) spec_2[c_1]']))
        self.assertEqual(set([p.serialize() for p in objs[SpeciesCoefficient].values()]),
                         set([p.serialize() for p in parts1 + parts2]))
        self.assertEqual(len(objs[Species]), 4)
        self.assertEqual(set(objs[Species].values()), set([p.species for p in parts1 + parts2]))

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

        # empty RHS
        parts, error = attr.deserialize('spec_2[c_1] ==>', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts]), set(
            ['(-1) spec_2[c_1]']))

        parts, error = attr.deserialize('[c_1]: spec_2 ==>', objs)
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

    def test_RateLawEquationAttribute_serialize(self):
        rxn = self.rxn_0
        rate_law = rxn.rate_laws[0]
        equation = rate_law.equation

        attr = RateLawEquationAttribute()
        self.assertEqual(attr.serialize(equation), equation.expression)

    def test_RateLawEquationAttribute_deserialize(self):
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
        }

        expression = 'k_cat * spec_0[c_0]'
        attr = RateLawEquationAttribute()
        equation1, error = attr.deserialize(expression, objs)
        self.assertEqual(error, None)
        self.assertEqual(equation1.expression, expression)
        self.assertEqual(equation1.modifiers[0].serialize(), 'spec_0[c_0]')
        self.assertEqual(list(objs[RateLawEquation].values()), [equation1])

    def test_objective_function_deserialize(self):

        objs = {
            Reaction: {
                'reaction_0': Reaction(id='reaction_0'),
                'reaction_1': Reaction(id='reaction_1'),
                'reaction_2': Reaction(id='reaction_2'),
            },
            BiomassReaction: {
                'biomass_reaction_0': BiomassReaction(id='biomass_reaction_0'),
                'biomass_reaction_1': BiomassReaction(id='biomass_reaction_1'),
            },
        }

        attr = ObjectiveFunction.Meta.attributes['expression']

        value = None
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, value, objs)
        self.assertTrue(of is None)
        self.assertTrue(invalid_attribute is None)

        value = "2*biomass_reaction_1 - pow( reaction_1, 2)"
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, value, objs)
        self.assertEqual(of.reactions[0], objs[Reaction]['reaction_1'])
        self.assertEqual(of.biomass_reactions[0], objs[BiomassReaction]['biomass_reaction_1'])

        objs[Reaction]['biomass_reaction_1'] = Reaction(id='biomass_reaction_1')
        value = "2*biomass_reaction_1 - pow( reaction_1, 2)"
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, value, objs)
        self.assertTrue(of is None)
        self.assertIn("id 'biomass_reaction_1' ambiguous between a Reaction and a BiomassReaction",
                      invalid_attribute.messages[0])

        del objs[Reaction]['biomass_reaction_1']
        value = "2*biomass_reaction_1 - pow( reaction_x, 2)"
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, value, objs)
        self.assertTrue(of is None)
        self.assertIn("id 'reaction_x' not a Reaction or a BiomassReaction identifier",
                      invalid_attribute.messages[0])

        value = "2*biomass_reaction_1 - pow( biomass_reaction_1, 2)"
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, value, objs)
        self.assertEqual(of.reactions, [])
        self.assertEqual(of.biomass_reactions[0], objs[BiomassReaction]['biomass_reaction_1'])
        self.assertEqual(len(of.biomass_reactions), 1)

    def test_objective_function_deserialize_invalid_ids(self):

        objs = {
            Reaction: {
                'pow': Reaction(id='reaction_0'),
                'reaction_x': Reaction(id='reaction_x'),
            },
            BiomassReaction: {
                'exp': BiomassReaction(id='biomass_reaction_0'),
            },
        }

        attr = ObjectiveFunction.Meta.attributes['expression']
        value = "2*exp - pow( reaction_x, 2)"
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, value, objs)
        self.assertTrue(of is None)
        self.assertIn("reaction id(s) {pow} ambiguous between a Reaction and a valid function",
                      invalid_attribute.messages[0])
        self.assertIn("reaction id(s) {exp} ambiguous between a BiomassReaction and a valid function",
                      invalid_attribute.messages[1])

    def test_objective_function_validate(self):

        objs = {
            Reaction: {
                'reaction_0': Reaction(id='reaction_0'),
                'reaction_1': Reaction(id='reaction_1'),
                'reaction_2': Reaction(id='reaction_2'),
            },
            BiomassReaction: {
                'biomass_reaction_0': BiomassReaction(id='biomass_reaction_0'),
                'biomass_reaction_1': BiomassReaction(id='biomass_reaction_1'),
            },
        }

        attr = ObjectiveFunction.Meta.attributes['expression']

        value = "2*biomass_reaction_1 - pow( reaction_1, 2)"
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, value, objs)
        self.assertTrue(invalid_attribute is None)
        rv = of.validate()
        self.assertTrue(rv is None)

        value = "2*biomass_reaction_1 - pow( reaction_1, 2"
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, value, objs)
        rv = of.validate()
        self.assertTrue(isinstance(rv, InvalidObject))
        self.assertEqual(rv.attributes[0].messages[0], "syntax error in expression '{}'".format(value))

        value = "2*biomass_reaction_1 - pow( 3*reaction_1, 2)"
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, value, objs)
        of.biomass_reactions = []
        rv = of.validate()
        self.assertEqual(rv.attributes[0].messages[0], "NameError in expression '{}'".format(value))

        obj_func = ObjectiveFunction(expression="'str' + 1.")
        self.assertNotEqual(obj_func.validate(), None)

    def test_validate(self):
        self.assertEqual(self.model.validate(), None)

    def test_sbml_data_exchange(self):
        # create an SBMLDocument that uses version 2 of the 'Flux Balance Constraints' extension
        document = create_sbml_doc_w_fbc()

        # Initialize the SBML document's model
        sbml_model = init_sbml_model(document)

        # Write a dFBA Submodel to an SBML document
        self.submdl_2.comments = 'test submodel comment'
        sbml_model = self.submdl_2.add_to_sbml_doc(document)
        self.assertEqual(sbml_model.getIdAttribute(), self.submdl_2.id)
        self.assertEqual(sbml_model.getName(), self.submdl_2.name)
        self.assertIn(self.submdl_2.comments, sbml_model.getNotesString())

        # Write Compartments to the SBML document
        self.comp_0.comments = 'test comment'
        sbml_compartment = self.comp_0.add_to_sbml_doc(document)
        self.assertTrue(sbml_compartment.hasRequiredAttributes())
        self.assertEqual(sbml_compartment.getIdAttribute(), self.comp_0.id)
        self.assertEqual(sbml_compartment.getName(), self.comp_0.name)
        self.assertEqual(sbml_compartment.getSize(), self.comp_0.initial_volume)
        self.assertIn(self.comp_0.comments, sbml_compartment.getNotesString())

        # Write species used by the submodel to the SBML document
        for species in self.submdl_2.get_species():
            sbml_species = species.add_to_sbml_doc(document)
            self.assertTrue(sbml_species.hasRequiredAttributes())
            self.assertEqual(sbml_species.getIdAttribute(), species.xml_id())
            self.assertEqual(sbml_species.getName(), species.species_type.name)
            self.assertEqual(sbml_species.getCompartment(), species.compartment.id)
            self.assertEqual(sbml_species.getInitialConcentration(), species.concentration.value)

        # Write reactions used by the submodel to an SBML document
        self.rxn_2.min_flux = 100
        self.rxn_2.max_flux = 200
        self.rxn_2.comments = 'comments'
        sbml_reaction = self.rxn_2.add_to_sbml_doc(document)
        self.assertTrue(sbml_reaction.hasRequiredAttributes())
        self.assertEqual(sbml_reaction.getIdAttribute(), self.rxn_2.id)
        self.assertEqual(sbml_reaction.getName(), self.rxn_2.name)
        self.assertEqual(sbml_reaction.getCompartment(), self.rxn_2.submodel.compartment.id)
        fbc_plugin = sbml_reaction.getPlugin('fbc')
        sbml_model = document.getModel()
        self.assertEqual(sbml_model.getParameter(fbc_plugin.getLowerFluxBound()).getValue(),
                         self.rxn_2.min_flux)
        self.assertEqual(sbml_model.getParameter(fbc_plugin.getUpperFluxBound()).getValue(),
                         self.rxn_2.max_flux)
        self.assertEqual(len(sbml_reaction.getListOfReactants()) + len(sbml_reaction.getListOfProducts()),
                         len(self.rxn_2.participants))
        for reactant in sbml_reaction.getListOfReactants():
            for participant in self.rxn_2.participants:
                if reactant.getSpecies() == participant.species.xml_id():
                    self.assertEqual(reactant.getStoichiometry(), -participant.coefficient)
        for product in sbml_reaction.getListOfProducts():
            for participant in self.rxn_2.participants:
                if product.getSpecies() == participant.species.xml_id():
                    self.assertEqual(product.getStoichiometry(), participant.coefficient)

        # Write the biomass reaction to the SBML document
        sbml_biomass_reaction = self.biomass_reaction.add_to_sbml_doc(document)
        self.assertTrue(sbml_biomass_reaction.hasRequiredAttributes())
        self.assertEqual(sbml_biomass_reaction.getIdAttribute(), self.biomass_reaction.id)
        self.assertEqual(sbml_biomass_reaction.getName(), self.biomass_reaction.name)
        self.assertIn(self.biomass_reaction.comments, sbml_biomass_reaction.getNotesString())
        fbc_plugin = sbml_biomass_reaction.getPlugin('fbc')
        sbml_model = document.getModel()
        self.assertEqual(sbml_model.getParameter(fbc_plugin.getLowerFluxBound()).getValue(), 0)
        self.assertEqual(sbml_model.getParameter(fbc_plugin.getUpperFluxBound()).getValue(),
                         float('inf'))
        self.assertEqual(len(sbml_biomass_reaction.getListOfReactants()) +
                         len(sbml_biomass_reaction.getListOfProducts()),
                         len(self.biomass_reaction.biomass_components))

        # Write parameters to the SBML document
        param = self.model.parameters.create(
            id='param_custom_units', name='param custom units',
            value=100, units='custom')
        param.submodel = self.model.submodels[0]
        self.parameters.append(param)

        for param in self.parameters:
            sbml_param = param.add_to_sbml_doc(document)
            self.assertTrue(sbml_param.hasRequiredAttributes())
            self.assertIn(param.id, sbml_param.getIdAttribute())
            self.assertEqual(sbml_param.getName(), param.name)
            self.assertEqual(sbml_param.getValue(), param.value)

        # Write an objective function to the model
        #   create objectiveFunction
        attr = ObjectiveFunction.Meta.attributes['expression']
        rxn_id = 'rxn_2'
        biomass_reaction_id = 'biomass_reaction_1'
        objs = {
            Reaction: {
                rxn_id: self.rxn_2,
            },
            BiomassReaction: {
                biomass_reaction_id: self.biomass_reaction},
        }
        (of, _) = ObjectiveFunction.deserialize(attr, 'biomass_reaction_1 + 2*rxn_2', objs)
        self.submdl_2.objective_function = of

        prepare_model = PrepareModel(self.model)
        (reactions, biomass_reactions) = prepare_model.parse_dfba_submodel_obj_func(self.submdl_2)
        PrepareModel.assign_linear_objective_fn(self.submdl_2, reactions, biomass_reactions)
        self.submdl_2.objective_function.linear = True

        #   write ObjectiveFunction to the model, and test
        sbml_objective = of.add_to_sbml_doc(document)
        self.assertEqual(wrap_libsbml(sbml_objective.getNumFluxObjectives, returns_int=True), 2)
        self.assertEqual(len(wrap_libsbml(sbml_objective.getListOfFluxObjectives)), 2)
        for flux_objective in wrap_libsbml(sbml_objective.getListOfFluxObjectives):
            if wrap_libsbml(flux_objective.getReaction) == rxn_id:
                self.assertEqual(wrap_libsbml(flux_objective.getCoefficient), 2.0)
            elif wrap_libsbml(flux_objective.getReaction) == biomass_reaction_id:
                self.assertEqual(wrap_libsbml(flux_objective.getCoefficient), 1.0)
            else:
                self.fail("reaction {} unexpected".format(wrap_libsbml(flux_objective.getReaction)))

        # Check the SBML document
        self.assertEqual(wrap_libsbml(get_SBML_compatibility_method(document)), 0)
        for i in range(document.checkConsistency()):
            print(document.getError(i).getShortMessage())
            print(document.getError(i).getMessage())
        self.assertEqual(wrap_libsbml(document.checkConsistency), 0)

        # exceptions
        obj_func = ObjectiveFunction(linear=False, submodels=[Submodel(id='Metabolism')])
        with pytest.warns(UserWarning):
            obj_func.add_to_sbml_doc(None)

    def test_objective_function_get_products(self):
        model = Model()
        submodel = model.submodels.create()
        species_type_0 = model.species_types.create(id='spec_0')
        species_type_1 = model.species_types.create(id='spec_1')
        species_type_2 = model.species_types.create(id='spec_2')
        compartment_0 = model.compartments.create(id='c_0')
        compartment_1 = model.compartments.create(id='c_1')
        compartment_2 = model.compartments.create(id='c_2')
        species_0 = Species(species_type=species_type_0, compartment=compartment_0)
        species_1 = Species(species_type=species_type_1, compartment=compartment_1)
        species_2 = Species(species_type=species_type_2, compartment=compartment_2)

        obj_func = submodel.objective_function = ObjectiveFunction(
            reactions=[
                Reaction(
                    reversible=True,
                    participants=[SpeciesCoefficient(species=species_0)],
                ),
            ],
            biomass_reactions=[
                BiomassReaction(
                    biomass_components=[BiomassComponent(coefficient=-1, species_type=species_type_1)],
                    compartment=compartment_1,
                ),
            ],
        )
        self.assertEqual(obj_func.get_products(), [species_0])

        obj_func = submodel.objective_function = ObjectiveFunction(
            reactions=[
                Reaction(
                    reversible=True,
                    participants=[SpeciesCoefficient(species=species_0)],
                ),
            ],
            biomass_reactions=[
                BiomassReaction(
                    biomass_components=[BiomassComponent(coefficient=-1, species_type=species_type_1)],
                    compartment=compartment_1,
                ),
                BiomassReaction(
                    biomass_components=[BiomassComponent(coefficient=1, species_type=species_type_2)],
                    compartment=compartment_2,
                ),
            ],
        )
        with self.assertRaisesRegexp(ValueError, 'does not belong to submodel'):
            obj_func.get_products()

    def test_function_validate(self):
        model = Model()
        model.species_types.create(id='A')
        model.species_types.create(id='BB')
        model.compartments.create(id='a')
        model.compartments.create(id='bb')
        model.observables.create(id='CCC')
        model.observables.create(id='DDD')

        func = model.functions.create(id='func', expression='CCC')
        self.assertEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='CCC + DDD')
        self.assertEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='CCC + 2 * DDD')
        self.assertEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='CCC + 2 * DDD > 3')
        self.assertEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='A[a] + BB[bb] + CCC > 3')
        self.assertNotEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='a[a] + BB[bb] + CCC > 3')
        self.assertNotEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='a[A] + BB[bb] + CCC > 3')
        self.assertNotEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='A[a] + BB[bb] + CC > 3')
        self.assertNotEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='"a"')
        self.assertNotEqual(func.validate(), None)

        func = model.functions.create(id='func', expression=' > 3')
        self.assertNotEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='_x')
        self.assertNotEqual(func.validate(), None)

        func = model.functions.create(id='func', expression='x() > 3')
        self.assertNotEqual(func.validate(), None)

    def test_stop_condition_validate(self):
        model = Model()
        model.species_types.create(id='A')
        model.species_types.create(id='BB')
        model.compartments.create(id='a')
        model.compartments.create(id='bb')
        model.observables.create(id='CCC')
        model.observables.create(id='DDD')

        cond = model.stop_conditions.create(id='cond', expression='CCC > 3')
        self.assertEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression='CCC + DDD > 3')
        self.assertEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression='CCC + 2 * DDD > 3')
        self.assertEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression='A[a] + BB[bb] + CCC > 3')
        self.assertNotEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression='a[a] + BB[bb] + CCC > 3')
        self.assertNotEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression='a[A] + BB[bb] + CCC > 3')
        self.assertNotEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression='A[a] + BB[bb] + CC > 3')
        self.assertNotEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression='CCC + 2 * DDD')
        self.assertNotEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression=' > 3')
        self.assertNotEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression='_x')
        self.assertNotEqual(cond.validate(), None)

        cond = model.stop_conditions.create(id='cond', expression='x() > 3')
        self.assertNotEqual(cond.validate(), None)


class TestCoreFromFile(unittest.TestCase):

    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_model.xlsx')

    def setUp(self):
        Submodel.objects.reset()
        # read and initialize a model
        self.model = Reader().run(self.MODEL_FILENAME)
        self.dfba_submodel = Submodel.objects.get_one(id='submodel_1')

    def test_get_ex_species(self):
        ex_compartment = self.model.compartments.get_one(id=EXTRACELLULAR_COMPARTMENT_ID)
        ex_species = self.dfba_submodel.get_species(compartment=ex_compartment)
        self.assertEqual(set(ex_species),
                         set(Species.get(['specie_1[e]', 'specie_2[e]'], self.dfba_submodel.get_species())))


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
