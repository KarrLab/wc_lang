""" Tests of core

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang.core import (Model, Taxon, TaxonRank, Submodel, Reaction, SpeciesType, SpeciesTypeType, Species, Compartment,
                          ReactionParticipant, Parameter, Reference, ReferenceType, CrossReference,
                          RateLaw, RateLawEquation, SubmodelAlgorithm, Concentration,
                          OneToOneSpeciesAttribute, ReactionParticipantsAttribute, RateLawEquationAttribute)
import unittest


class TestCore(unittest.TestCase):

    def setUp(self):
        self.model = mdl = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1b')

        mdl.taxon = Taxon(id='taxon', name='test taxon', rank=TaxonRank.species)

        self.comp_0 = comp_0 = mdl.compartments.create(id='comp_0', name='compartment 0', initial_volume=1.25)
        self.comp_1 = comp_1 = mdl.compartments.create(id='comp_1', name='compartment 1', initial_volume=2.5)
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

        self.submdl_0 = submdl_0 = mdl.submodels.create(
            id='submodel_0', name='submodel 0', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_1 = submdl_1 = mdl.submodels.create(
            id='submodel_1', name='submodel 1', algorithm=SubmodelAlgorithm.ssa)
        self.submdl_2 = submdl_2 = mdl.submodels.create(
            id='submodel_2', name='submodel 2', algorithm=SubmodelAlgorithm.dfba)
        self.submodels = submodels = [submdl_0, submdl_1, submdl_2]

        self.rxn_0 = rxn_0 = submdl_0.reactions.create(id='rxn_0', name='reaction 0')
        rxn_0.participants.create(species=species[0], coefficient=-2)
        rxn_0.participants.create(species=species[1], coefficient=-3.5)
        rxn_0.participants.create(species=species[2], coefficient=1)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[5].serialize()),
            modifiers=set(species[5:6]))
        rate_law_0 = rxn_0.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_1 = rxn_1 = submdl_1.reactions.create(id='rxn_1', name='reaction 1')
        rxn_1.participants.create(species=species[0], coefficient=-2)
        rxn_1.participants.create(species=species[1], coefficient=-3)
        rxn_1.participants.create(species=species[3], coefficient=2)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[6].serialize()),
            modifiers=set(species[6:7]))
        rate_law_1 = rxn_1.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_2 = rxn_2 = submdl_2.reactions.create(id='rxn_2', name='reaction 2')
        rxn_2.participants.create(species=species[0], coefficient=-2)
        rxn_2.participants.create(species=species[1], coefficient=-3)
        rxn_2.participants.create(species=species[4], coefficient=1)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[7].serialize()),
            modifiers=set(species[7:8]))
        rate_law_2 = rxn_2.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.reactions = [rxn_0, rxn_1, rxn_2]
        self.rate_laws = [rate_law_0, rate_law_1, rate_law_2]

        self.parameters = parameters = []
        self.references = references = []
        self.cross_references = cross_references = []
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

            x_ref = ref.cross_references.create(database='x', id='y' * (i + 1),
                                                url='http://x.com/{}'.format('y' * (i + 1)))
            cross_references.append(x_ref)

    def test_reverse_references(self):
        mdl = self.model

        self.assertEqual(mdl.submodels, set(self.submodels))
        self.assertEqual(mdl.compartments, set(self.compartments))
        self.assertEqual(mdl.species_types, set(self.species_types))
        self.assertEqual(mdl.parameters, set(self.parameters))

        # submodel
        for reaction, submodel in zip(self.reactions, self.submodels):
            self.assertEqual(submodel.reactions, set([reaction]))

        for submodel in self.submodels:
            self.assertEqual(submodel.cross_references, set())
            self.assertEqual(submodel.references, set())

        # compartment
        self.assertEqual(self.compartments[0].species, set(self.species[0:3] + self.species[4:]))
        self.assertEqual(self.compartments[1].species, set([self.species[3]]))

        for compartment in self.compartments:
            self.assertEqual(compartment.cross_references, set())
            self.assertEqual(compartment.references, set())

        # species type
        for species_type, species in zip(self.species_types, self.species):
            self.assertEqual(species_type.species, set([species]))

        for species_type in self.species_types:
            self.assertEqual(species_type.cross_references, set())
            self.assertEqual(species_type.references, set())

        # specie
        for species_type, species in zip(self.species_types, self.species):
            self.assertEqual(species.species_type, species_type)

        for i in range(len(self.species)):
            if i != 3:
                self.assertEqual(self.species[i].compartment, self.compartments[0])
            else:
                self.assertEqual(self.species[i].compartment, self.compartments[1])

        self.assertEqual(len(self.species[0].reaction_participants), 3)
        self.assertEqual(len(self.species[1].reaction_participants), 3)
        self.assertEqual(len(self.species[2].reaction_participants), 1)
        self.assertEqual(len(self.species[3].reaction_participants), 1)
        self.assertEqual(len(self.species[4].reaction_participants), 1)
        self.assertEqual(len(self.species[5].reaction_participants), 0)
        self.assertEqual(len(self.species[6].reaction_participants), 0)
        self.assertEqual(len(self.species[7].reaction_participants), 0)

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

        self.assertEqual(list(self.reactions[0].rate_laws)[0].equation.modifiers, set(self.species[5:6]))
        self.assertEqual(list(self.reactions[1].rate_laws)[0].equation.modifiers, set(self.species[6:7]))
        self.assertEqual(list(self.reactions[2].rate_laws)[0].equation.modifiers, set(self.species[7:8]))

        for reaction in self.reactions:
            self.assertEqual(reaction.cross_references, set())
            self.assertEqual(reaction.references, set())
            self.assertEqual(len(reaction.rate_laws), 1)

        # parameters
        for reference, parameter in zip(self.references, self.parameters):
            self.assertEqual(parameter.references, set([reference]))

        for parameter in self.parameters:
            self.assertEqual(parameter.model, mdl)

        # references
        for reference, parameter in zip(self.references, self.parameters):
            self.assertEqual(reference.parameters, set([parameter]))
            self.assertEqual(parameter.references, set([reference]))

        for reference, cross_reference in zip(self.references, self.cross_references):
            self.assertEqual(reference.cross_references, set([cross_reference]))
            self.assertEqual(cross_reference.reference, reference)

        # reaction participant
        for species in self.species[0:5]:
            self.assertEqual(set(x.species for x in species.reaction_participants), set([species]))

        for reaction in self.reactions:
            for part in reaction.participants:
                self.assertIn(reaction, part.reactions)
            self.assertEqual(set(x.reaction for x in reaction.rate_laws), set([reaction]))

        # cross references
        for reference, cross_reference in zip(self.references, self.cross_references):
            self.assertEqual(reference.cross_references, set([cross_reference]))
            self.assertEqual(cross_reference.reference, reference)

    def test_taxon_rank_class(self):
        self.assertEqual(TaxonRank['class'], TaxonRank['classis'])

    def test_model_get_species(self):
        self.assertEqual(self.model.get_species(), set(self.species))

    def test_submodel_get_species(self):
        species = self.species
        self.assertEqual(self.submdl_0.get_species(), set([
            species[0], species[1], species[2], species[5],
        ]))
        self.assertEqual(self.submdl_1.get_species(), set([
            species[0], species[1], species[3], species[6],
        ]))
        self.assertEqual(self.submdl_2.get_species(), set([
            species[0], species[1], species[4], species[7],
        ]))

    def test_reaction_get_species(self):
        species = self.species
        self.assertEqual(self.rxn_0.get_species(), set([
            species[0], species[1], species[2], species[5],
        ]))
        self.assertEqual(self.rxn_1.get_species(), set([
            species[0], species[1], species[3], species[6],
        ]))
        self.assertEqual(self.rxn_2.get_species(), set([
            species[0], species[1], species[4], species[7],
        ]))

    def test_get_components(self):
        mdl = self.model

        self.assertEqual(mdl.get_compartments(), set(self.compartments))
        self.assertEqual(mdl.get_species_types(), set(self.species_types))
        self.assertEqual(mdl.get_submodels(), set(self.submodels))
        self.assertEqual(mdl.get_species(), set(self.species))
        self.assertEqual(mdl.get_concentrations(), set(self.concentrations))
        self.assertEqual(mdl.get_reactions(), set(self.reactions))
        self.assertEqual(mdl.get_rate_laws(), set(self.rate_laws))
        self.assertEqual(mdl.get_parameters(), set(self.parameters))
        self.assertEqual(mdl.get_references(), set(self.references))

    def test_get_component(self):
        model = self.model

        self.assertEqual(model.get_component('compartment', 'comp_0'), self.comp_0)
        self.assertEqual(model.get_component('species_type', 'spec_type_1'), self.species_types[1])
        self.assertEqual(model.get_component('submodel', 'submodel_1'), self.submdl_1)
        self.assertEqual(model.get_component('reaction', 'rxn_1'), self.rxn_1)
        self.assertEqual(model.get_component('parameter', 'param_2'), self.parameters[2])
        self.assertEqual(model.get_component('reference', 'ref_1'), self.references[1])

        self.assertEqual(model.get_component('reaction', 'rxn_3'), None)

    def test_species_type_is_carbon_containing(self):
        self.assertFalse(self.species_types[0].is_carbon_containing())
        self.assertTrue(self.species_types[1].is_carbon_containing())

    def test_species_serialize(self):
        self.assertEqual(self.species[0].serialize(), 'spec_type_0[comp_0]')
        self.assertEqual(self.species[1].serialize(), 'spec_type_1[comp_0]')
        self.assertEqual(self.species[2].serialize(), 'spec_type_2[comp_0]')
        self.assertEqual(self.species[3].serialize(), 'spec_type_3[comp_1]')

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

    def test_concentration_serialize(self):
        self.assertEqual(self.concentrations[0].serialize(), 'spec_type_0[comp_0]')
        self.assertEqual(self.concentrations[1].serialize(), 'spec_type_1[comp_0]')
        self.assertEqual(self.concentrations[2].serialize(), 'spec_type_2[comp_0]')
        self.assertEqual(self.concentrations[3].serialize(), 'spec_type_3[comp_1]')

    def test_reaction_participant_serialize(self):
        self.assertEqual(set([part.serialize() for part in self.rxn_0.participants]), set([
            '(2) spec_type_0[comp_0]', '(3.5) spec_type_1[comp_0]', 'spec_type_2[comp_0]'
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
        part0, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part0.coefficient, 1)
        self.assertEqual(part0.species.serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[ReactionParticipant].values()), set([part0]))
        self.assertEqual(set(objs[Species].values()), set([part0.species]))

        val = '(2) spec_0[c_0]'
        part1, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part1.coefficient, 2)
        self.assertEqual(part1.species.serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[ReactionParticipant].values()), set([part0, part1]))
        self.assertEqual(set(objs[Species].values()), set([part0.species, part1.species]))

        val = '(2.) spec_0[c_1]'
        part2, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part2.coefficient, 2)
        self.assertEqual(part2.species.serialize(), 'spec_0[c_1]')
        self.assertEqual(set(objs[ReactionParticipant].values()), set([part0, part1, part2]))
        self.assertEqual(set(objs[Species].values()), set([part0.species, part1.species, part2.species]))

        val = '(2.5) spec_0[c_0]'
        part3, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part3.coefficient, 2.5)
        self.assertEqual(part3.species.serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[ReactionParticipant].values()), set([part0, part1, part2, part3]))
        self.assertEqual(set(objs[Species].values()), set([part0.species, part1.species, part2.species, part3.species]))

        val = '(.5) spec_0[c_0]'
        part4, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertEqual(error, None)
        self.assertEqual(part4.coefficient, 0.5)
        self.assertEqual(part4.species.serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[ReactionParticipant].values()), set([part0, part1, part2, part3, part4]))
        self.assertEqual(set(objs[Species].values()), set(
            [part0.species, part1.species, part2.species, part3.species, part4.species]))

        val = '(1) spec_1'
        part5, error = ReactionParticipant.deserialize(attr, val, objs, compartment=objs[Compartment]['c_0'])
        self.assertEqual(error, None)
        self.assertEqual(part5.coefficient, 1)
        self.assertEqual(part5.species.serialize(), 'spec_1[c_0]')
        self.assertEqual(set(objs[ReactionParticipant].values()), set([part0, part1, part2, part3, part4, part5]))
        self.assertEqual(set(objs[Species].values()), set(
            [part0.species, part1.species, part2.species, part3.species, part4.species, part5.species]))

        # negative examples
        val = '(-1) spec_0[c_0]'
        part6, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        val = '(1) spec_0'
        part6, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        val = '(1.1.) spec_0[c_0]'
        part6, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        val = ' spec_0[c_0]'
        part6, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        val = ' spec_3[c_0]'
        part6, error = ReactionParticipant.deserialize(attr, val, objs)
        self.assertNotEqual(error, None)
        self.assertEqual(part6, None)

        self.assertEqual(set(objs[ReactionParticipant].values()), set([part0, part1, part2, part3, part4, part5]))
        self.assertEqual(set(objs[Species].values()), set(
            [part0.species, part1.species, part2.species, part3.species, part4.species, part5.species]))

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
        self.assertEqual(list(equation1.modifiers)[0].serialize(), 'spec_0[c_0]')
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

        expression = 'k_cat * spec_0[c_0]'
        equation = RateLawEquation(
            expression=expression,
            modifiers=[
                Species(species_type=species_types[0], compartment=compartments[0])
            ])
        self.assertNotEqual(equation.validate(), None)

        expression = 'k_cat * spec_0[c_0]'
        equation = RateLawEquation(
            rate_law=RateLaw(k_cat=1, k_m=1),
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
        self.assertEqual(list(self.rxn_0.rate_laws)[0].equation.modifiers, set([self.species[5]]))
        self.assertEqual(list(self.rxn_1.rate_laws)[0].equation.modifiers, set([self.species[6]]))
        self.assertEqual(list(self.rxn_2.rate_laws)[0].equation.modifiers, set([self.species[7]]))

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
            Parameter(id='a', submodels=[Submodel()]),
            Parameter(id='a', submodels=[Submodel()]),
        ]
        self.assertEqual(Parameter.validate_unique(params), None)

    def test_cross_reference_serialize(self):
        self.assertEqual(self.cross_references[0].serialize(), '{}: {}'.format('x', 'y'))
        self.assertEqual(self.cross_references[1].serialize(), '{}: {}'.format('x', 'yy'))
        self.assertEqual(self.cross_references[2].serialize(), '{}: {}'.format('x', 'yyy'))

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
        self.assertEqual(set(objs[Species].values()), set([species0]))

    def test_ReactionParticipantsAttribute_serialize(self):
        attr = ReactionParticipantsAttribute()
        self.assertEqual(attr.serialize(self.rxn_0.participants),
                         '[comp_0]: (2) spec_type_0 + (3.5) spec_type_1 ==> spec_type_2')
        self.assertEqual(attr.serialize(self.rxn_1.participants),
                         '(2) spec_type_0[comp_0] + (3) spec_type_1[comp_0] ==> (2) spec_type_3[comp_1]')

    def test_ReactionParticipantsAttribute_deserialize(self):
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

        attr = ReactionParticipantsAttribute()

        parts1, error = attr.deserialize('[c_0]: (2) spec_0 + (3.5) spec_1 ==> spec_2', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts1]), set(
            ['(2) spec_0[c_0]', '(3.5) spec_1[c_0]', 'spec_2[c_0]']))
        self.assertEqual(len(objs[ReactionParticipant]), 3)
        self.assertEqual(set(objs[ReactionParticipant].values()), parts1)
        self.assertEqual(len(objs[Species]), 3)
        self.assertEqual(set(objs[Species].values()), set([p.species for p in parts1]))

        parts2, error = attr.deserialize(
            '(2) spec_0[c_0] + (3) spec_1[c_0] ==> (2) spec_2[c_1]', objs)
        self.assertEqual(error, None)
        self.assertEqual(set([p.serialize() for p in parts2]), set(
            ['(2) spec_0[c_0]', '(3) spec_1[c_0]', '(2) spec_2[c_1]']))
        self.assertEqual(set([p.serialize() for p in objs[ReactionParticipant].values()]),
                         set([p.serialize() for p in parts1 | parts2]))
        self.assertEqual(len(objs[Species]), 4)
        self.assertEqual(set(objs[Species].values()), set([p.species for p in parts1 | parts2]))

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

    def test_RateLawEquationAttribute_serialize(self):
        rxn = self.rxn_0
        rate_law = list(rxn.rate_laws)[0]
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
        self.assertEqual(list(equation1.modifiers)[0].serialize(), 'spec_0[c_0]')
        self.assertEqual(set(objs[RateLawEquation].values()), set([equation1]))

    def test_validate(self):
        self.assertEqual(self.model.validate(), None)
