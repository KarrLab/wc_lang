""" Tests of merging

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-08
:Copyright: 2019, Karr Lab
:License: MIT
"""

import unittest
import wc_lang.util
from wc_lang.core import (Model, Taxon, TaxonRank, Environment,
                          Submodel, Reaction, SpeciesCoefficient,
                          ObservableExpression, FunctionExpression,
                          RateLawExpression, DfbaObjectiveExpression,
                          StopConditionExpression,
                          Evidence, Identifier, Reference)
from wc_onto import onto
from wc_utils.util.units import unit_registry


class MergeTestCase(unittest.TestCase):

    def gen_models(self):
        # model
        models = [
            Model(id='model_1', name='model 1', version='0.0.1', wc_lang_version='0.0.1'),
            Model(id='model_2', name='model 2', version='0.0.2', wc_lang_version='0.0.2'),
            Model(id='model_1', name='model 1', version='0.0.1', wc_lang_version='0.0.1'),  # left merge
            Model(id='model_2', name='model 2', version='0.0.2', wc_lang_version='0.0.2'),  # right merge
        ]

        # taxon
        models[0].taxon = Taxon(id='taxon', name='taxon 1', rank=TaxonRank.species, comments='comments 1')
        models[1].taxon = Taxon(id='taxon', name='taxon 1', rank=TaxonRank.species, comments='comments 1')
        models[2].taxon = Taxon(id='taxon', name='taxon 1', rank=TaxonRank.species, comments='comments 1')
        models[3].taxon = Taxon(id='taxon', name='taxon 1', rank=TaxonRank.species, comments='comments 1')

        # environment
        models[0].env = Environment(id='env', comments='comments 1')
        models[1].env = Environment(id='env', comments='comments 2')
        models[2].env = Environment(id='env', comments='comments 1\n\ncomments 2')
        models[3].env = Environment(id='env', comments='comments 2\n\ncomments 1')

        # compartments
        models[0].compartments.create(id='c_0', mean_init_volume=1.)
        models[0].compartments.create(id='c_1', mean_init_volume=2., parent_compartment=models[0].compartments[0])

        models[1].compartments.create(id='c_0', mean_init_volume=1.)
        models[1].compartments.create(id='c_2', mean_init_volume=3.)

        models[2].compartments.create(id='c_0', mean_init_volume=1.)
        models[2].compartments.create(id='c_1', mean_init_volume=2., parent_compartment=models[2].compartments[0])
        models[2].compartments.create(id='c_2', mean_init_volume=3.)

        models[3].compartments.create(id='c_0', mean_init_volume=1.)
        models[3].compartments.create(id='c_1', mean_init_volume=2., parent_compartment=models[3].compartments[0])
        models[3].compartments.create(id='c_2', mean_init_volume=3.)

        for model in models:
            for comp in model.compartments:
                comp.geometry = None

        # species types
        models[0].species_types.create(id='st_0')
        models[0].species_types.create(id='st_1')
        models[0].species_types.create(id='st_2')

        models[1].species_types.create(id='st_0')
        models[1].species_types.create(id='st_1')
        models[1].species_types.create(id='st_3')

        models[2].species_types.create(id='st_0')
        models[2].species_types.create(id='st_1')
        models[2].species_types.create(id='st_2')
        models[2].species_types.create(id='st_3')

        models[3].species_types.create(id='st_0')
        models[3].species_types.create(id='st_1')
        models[3].species_types.create(id='st_2')
        models[3].species_types.create(id='st_3')

        for model in models:
            for species_type in model.species_types:
                species_type.charge = 0

        # species
        models[0].species.create(species_type=models[0].species_types[0], compartment=models[0].compartments[0])
        models[0].species.create(species_type=models[0].species_types[0], compartment=models[0].compartments[1])
        models[0].species.create(species_type=models[0].species_types[2], compartment=models[0].compartments[1])

        models[1].species.create(species_type=models[1].species_types[0], compartment=models[1].compartments[0])
        models[1].species.create(species_type=models[1].species_types[1], compartment=models[1].compartments[0])
        models[1].species.create(species_type=models[1].species_types[2], compartment=models[1].compartments[1])

        models[2].species.create(species_type=models[2].species_types[0], compartment=models[2].compartments[0])
        models[2].species.create(species_type=models[2].species_types[0], compartment=models[2].compartments[1])
        models[2].species.create(species_type=models[2].species_types[2], compartment=models[2].compartments[1])
        models[2].species.create(species_type=models[2].species_types[1], compartment=models[2].compartments[0])
        models[2].species.create(species_type=models[2].species_types[3], compartment=models[2].compartments[2])

        models[3].species.create(species_type=models[3].species_types[0], compartment=models[3].compartments[0])
        models[3].species.create(species_type=models[3].species_types[0], compartment=models[3].compartments[1])
        models[3].species.create(species_type=models[3].species_types[2], compartment=models[3].compartments[1])
        models[3].species.create(species_type=models[3].species_types[1], compartment=models[3].compartments[0])
        models[3].species.create(species_type=models[3].species_types[3], compartment=models[3].compartments[2])

        for model in models:
            for species in model.species:
                species.id = species.gen_id()

        # distribution initial conditions
        models[0].distribution_init_concentrations.create(species=models[0].species[0])
        models[0].distribution_init_concentrations.create(species=models[0].species[1])

        models[1].distribution_init_concentrations.create(species=models[1].species[0])
        models[1].distribution_init_concentrations.create(species=models[1].species[2])

        models[2].distribution_init_concentrations.create(species=models[2].species[0])
        models[2].distribution_init_concentrations.create(species=models[2].species[1])
        models[2].distribution_init_concentrations.create(species=models[2].species[4])

        models[3].distribution_init_concentrations.create(species=models[3].species[0])
        models[3].distribution_init_concentrations.create(species=models[3].species[1])
        models[3].distribution_init_concentrations.create(species=models[3].species[4])

        for model in models:
            for obj in model.distribution_init_concentrations:
                obj.id = obj.gen_id()

        # parameters
        models[0].parameters.create(id='p_0', units=unit_registry.parse_units('s^-1'))
        models[0].parameters.create(id='p_1', density_compartment=models[0].compartments[0], units=unit_registry.parse_units('g l^-1'))
        models[0].parameters.create(id='p_2', units=unit_registry.parse_units('g l^-1'))
        models[0].parameters.create(id='p_3', units=unit_registry.parse_units('s^-1'))
        models[0].parameters.create(id='p_4', density_compartment=models[0].compartments[1], units=unit_registry.parse_units('g l^-1'))

        models[1].parameters.create(id='p_0', units=unit_registry.parse_units('s^-1'))
        models[1].parameters.create(id='p_1', units=unit_registry.parse_units('g l^-1'))
        models[1].parameters.create(id='p_2', density_compartment=models[1].compartments[1], units=unit_registry.parse_units('g l^-1'))
        models[1].parameters.create(id='p_4', units=unit_registry.parse_units('g l^-1'))
        models[1].parameters.create(id='p_5', units=unit_registry.parse_units('s^-1'))
        models[1].parameters.create(id='p_6', units=unit_registry.parse_units('s^-1'))

        models[2].parameters.create(id='p_0', units=unit_registry.parse_units('s^-1'))
        models[2].parameters.create(id='p_1', density_compartment=models[2].compartments[0], units=unit_registry.parse_units('g l^-1'))
        models[2].parameters.create(id='p_2', density_compartment=models[2].compartments[2], units=unit_registry.parse_units('g l^-1'))
        models[2].parameters.create(id='p_3', units=unit_registry.parse_units('s^-1'))
        models[2].parameters.create(id='p_4', density_compartment=models[2].compartments[1], units=unit_registry.parse_units('g l^-1'))
        models[2].parameters.create(id='p_5', units=unit_registry.parse_units('s^-1'))
        models[2].parameters.create(id='p_6', units=unit_registry.parse_units('s^-1'))

        models[3].parameters.create(id='p_0', units=unit_registry.parse_units('s^-1'))
        models[3].parameters.create(id='p_1', density_compartment=models[3].compartments[0], units=unit_registry.parse_units('g l^-1'))
        models[3].parameters.create(id='p_2', density_compartment=models[3].compartments[2], units=unit_registry.parse_units('g l^-1'))
        models[3].parameters.create(id='p_3', units=unit_registry.parse_units('s^-1'))
        models[3].parameters.create(id='p_4', density_compartment=models[3].compartments[1], units=unit_registry.parse_units('g l^-1'))
        models[3].parameters.create(id='p_5', units=unit_registry.parse_units('s^-1'))
        models[3].parameters.create(id='p_6', units=unit_registry.parse_units('s^-1'))

        # submodels
        models[0].submodels.create(id='submodel_11', framework=onto['WC:stochastic_simulation_algorithm'])
        models[0].submodels.create(id='submodel_12', framework=onto['WC:dynamic_flux_balance_analysis'])

        models[1].submodels.create(id='submodel_21', framework=onto['WC:stochastic_simulation_algorithm'])
        models[1].submodels.create(id='submodel_22', framework=onto['WC:dynamic_flux_balance_analysis'])

        models[2].submodels.create(id='submodel_11', framework=onto['WC:stochastic_simulation_algorithm'])
        models[2].submodels.create(id='submodel_12', framework=onto['WC:dynamic_flux_balance_analysis'])
        models[2].submodels.create(id='submodel_21', framework=onto['WC:stochastic_simulation_algorithm'])
        models[2].submodels.create(id='submodel_22', framework=onto['WC:dynamic_flux_balance_analysis'])

        models[3].submodels.create(id='submodel_11', framework=onto['WC:stochastic_simulation_algorithm'])
        models[3].submodels.create(id='submodel_12', framework=onto['WC:dynamic_flux_balance_analysis'])
        models[3].submodels.create(id='submodel_21', framework=onto['WC:stochastic_simulation_algorithm'])
        models[3].submodels.create(id='submodel_22', framework=onto['WC:dynamic_flux_balance_analysis'])

        # reactions
        models[0].reactions.create(id='rxn_11', submodel=models[0].submodels[0])
        models[0].reactions.create(id='rxn_12', submodel=models[0].submodels[0])

        models[1].reactions.create(id='rxn_31', submodel=models[1].submodels[0])
        models[1].reactions.create(id='rxn_41', submodel=models[1].submodels[1])

        models[2].reactions.create(id='rxn_11', submodel=models[2].submodels[0])
        models[2].reactions.create(id='rxn_12', submodel=models[2].submodels[0])
        models[2].reactions.create(id='rxn_31', submodel=models[2].submodels[2])
        models[2].reactions.create(id='rxn_41', submodel=models[2].submodels[3])

        models[3].reactions.create(id='rxn_11', submodel=models[3].submodels[0])
        models[3].reactions.create(id='rxn_12', submodel=models[3].submodels[0])
        models[3].reactions.create(id='rxn_31', submodel=models[3].submodels[2])
        models[3].reactions.create(id='rxn_41', submodel=models[3].submodels[3])

        # Reaction participants
        sc_0 = SpeciesCoefficient(species=models[0].species[0], coefficient=1.)
        sc_1 = SpeciesCoefficient(species=models[0].species[0], coefficient=2.)
        sc_2 = SpeciesCoefficient(species=models[0].species[1], coefficient=1.)
        sc_3 = SpeciesCoefficient(species=models[0].species[1], coefficient=2.)
        models[0].reactions[0].participants.append(sc_0)
        models[0].reactions[0].participants.append(sc_1)
        models[0].reactions[0].participants.append(sc_2)
        models[0].reactions[1].participants.append(sc_1)
        models[0].reactions[1].participants.append(sc_2)
        models[0].reactions[1].participants.append(sc_3)

        sc_0 = SpeciesCoefficient(species=models[1].species[0], coefficient=1.)
        sc_1 = SpeciesCoefficient(species=models[1].species[0], coefficient=3.)
        sc_2 = SpeciesCoefficient(species=models[1].species[1], coefficient=1.)
        sc_3 = SpeciesCoefficient(species=models[1].species[1], coefficient=3.)
        models[1].reactions[0].participants.append(sc_0)
        models[1].reactions[0].participants.append(sc_1)
        models[1].reactions[0].participants.append(sc_2)
        models[1].reactions[1].participants.append(sc_1)
        models[1].reactions[1].participants.append(sc_2)
        models[1].reactions[1].participants.append(sc_3)

        sc_0 = SpeciesCoefficient(species=models[2].species[0], coefficient=1.)
        sc_1 = SpeciesCoefficient(species=models[2].species[0], coefficient=2.)
        sc_2 = SpeciesCoefficient(species=models[2].species[1], coefficient=1.)
        sc_3 = SpeciesCoefficient(species=models[2].species[1], coefficient=2.)
        sc_4 = SpeciesCoefficient(species=models[2].species[0], coefficient=3.)
        sc_5 = SpeciesCoefficient(species=models[2].species[3], coefficient=1.)
        sc_6 = SpeciesCoefficient(species=models[2].species[3], coefficient=3.)
        models[2].reactions[0].participants.append(sc_0)
        models[2].reactions[0].participants.append(sc_1)
        models[2].reactions[0].participants.append(sc_2)
        models[2].reactions[1].participants.append(sc_1)
        models[2].reactions[1].participants.append(sc_2)
        models[2].reactions[1].participants.append(sc_3)
        models[2].reactions[2].participants.append(sc_0)
        models[2].reactions[2].participants.append(sc_4)
        models[2].reactions[2].participants.append(sc_5)
        models[2].reactions[3].participants.append(sc_4)
        models[2].reactions[3].participants.append(sc_5)
        models[2].reactions[3].participants.append(sc_6)

        sc_0 = SpeciesCoefficient(species=models[3].species[0], coefficient=1.)
        sc_1 = SpeciesCoefficient(species=models[3].species[0], coefficient=2.)
        sc_2 = SpeciesCoefficient(species=models[3].species[1], coefficient=1.)
        sc_3 = SpeciesCoefficient(species=models[3].species[1], coefficient=2.)
        sc_4 = SpeciesCoefficient(species=models[3].species[0], coefficient=3.)
        sc_5 = SpeciesCoefficient(species=models[3].species[3], coefficient=1.)
        sc_6 = SpeciesCoefficient(species=models[3].species[3], coefficient=3.)
        models[3].reactions[0].participants.append(sc_0)
        models[3].reactions[0].participants.append(sc_1)
        models[3].reactions[0].participants.append(sc_2)
        models[3].reactions[1].participants.append(sc_1)
        models[3].reactions[1].participants.append(sc_2)
        models[3].reactions[1].participants.append(sc_3)
        models[3].reactions[2].participants.append(sc_0)
        models[3].reactions[2].participants.append(sc_4)
        models[3].reactions[2].participants.append(sc_5)
        models[3].reactions[3].participants.append(sc_4)
        models[3].reactions[3].participants.append(sc_5)
        models[3].reactions[3].participants.append(sc_6)

        # rate laws
        models[0].rate_laws.create(reaction=models[0].reactions[0])
        models[0].rate_laws.create(reaction=models[0].reactions[1])

        models[1].rate_laws.create(reaction=models[1].reactions[0])
        models[1].rate_laws.create(reaction=models[1].reactions[1])

        models[2].rate_laws.create(reaction=models[2].reactions[0])
        models[2].rate_laws.create(reaction=models[2].reactions[1])
        models[2].rate_laws.create(reaction=models[2].reactions[2])
        models[2].rate_laws.create(reaction=models[2].reactions[3])

        models[3].rate_laws.create(reaction=models[3].reactions[0])
        models[3].rate_laws.create(reaction=models[3].reactions[1])
        models[3].rate_laws.create(reaction=models[3].reactions[2])
        models[3].rate_laws.create(reaction=models[3].reactions[3])

        for model in models:
            for rl in model.rate_laws:
                rl.id = rl.gen_id()

        objs = models[0].gen_serialized_val_obj_map()
        models[0].rate_laws[0].expression, _ = RateLawExpression.deserialize(
            '{}'.format(models[0].parameters[0].id), objs)
        models[0].rate_laws[1].expression, _ = RateLawExpression.deserialize(
            '2 * {}'.format(models[0].parameters[0].id), objs)

        objs = models[1].gen_serialized_val_obj_map()
        models[1].rate_laws[0].expression, _ = RateLawExpression.deserialize(
            '{}'.format(models[1].parameters[0].id), objs)
        models[1].rate_laws[1].expression, _ = RateLawExpression.deserialize(
            '2 * {}'.format(models[1].parameters[0].id), objs)

        objs = models[2].gen_serialized_val_obj_map()
        models[2].rate_laws[0].expression, _ = RateLawExpression.deserialize(
            '{}'.format(models[2].parameters[0].id), objs)
        models[2].rate_laws[1].expression, _ = RateLawExpression.deserialize(
            '2 * {}'.format(models[2].parameters[0].id), objs)
        models[2].rate_laws[2].expression, _ = RateLawExpression.deserialize(
            '{}'.format(models[2].parameters[0].id), objs)
        models[2].rate_laws[3].expression, _ = RateLawExpression.deserialize(
            '2 * {}'.format(models[2].parameters[0].id), objs)

        objs = models[3].gen_serialized_val_obj_map()
        models[3].rate_laws[0].expression, _ = RateLawExpression.deserialize(
            '{}'.format(models[3].parameters[0].id), objs)
        models[3].rate_laws[1].expression, _ = RateLawExpression.deserialize(
            '2 * {}'.format(models[3].parameters[0].id), objs)
        models[3].rate_laws[2].expression, _ = RateLawExpression.deserialize(
            '{}'.format(models[3].parameters[0].id), objs)
        models[3].rate_laws[3].expression, _ = RateLawExpression.deserialize(
            '2 * {}'.format(models[3].parameters[0].id), objs)

        # observables
        models[0].observables.create(id='obs_1')
        models[0].observables.create(id='obs_2')
        models[0].observables.create(id='obs_3')
        models[0].observables.create(id='obs_4')

        models[1].observables.create(id='obs_1')
        models[1].observables.create(id='obs_5')
        models[1].observables.create(id='obs_6')
        models[1].observables.create(id='obs_7')

        models[2].observables.create(id='obs_1')
        models[2].observables.create(id='obs_2')
        models[2].observables.create(id='obs_3')
        models[2].observables.create(id='obs_4')
        models[2].observables.create(id='obs_5')
        models[2].observables.create(id='obs_6')
        models[2].observables.create(id='obs_7')

        models[3].observables.create(id='obs_1')
        models[3].observables.create(id='obs_2')
        models[3].observables.create(id='obs_3')
        models[3].observables.create(id='obs_4')
        models[3].observables.create(id='obs_5')
        models[3].observables.create(id='obs_6')
        models[3].observables.create(id='obs_7')

        objs = models[0].gen_serialized_val_obj_map()
        models[0].observables[0].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[0].species[0].id), objs)
        models[0].observables[1].expression, _ = ObservableExpression.deserialize('2 * {}'.format(models[0].species[0].id), objs)
        models[0].observables[2].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[0].species[1].id), objs)
        models[0].observables[3].expression, _ = ObservableExpression.deserialize('2 * {}'.format(models[0].species[1].id), objs)

        objs = models[1].gen_serialized_val_obj_map()
        models[1].observables[0].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[1].species[0].id), objs)
        models[1].observables[1].expression, _ = ObservableExpression.deserialize('3 * {}'.format(models[1].species[0].id), objs)
        models[1].observables[2].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[1].species[1].id), objs)
        models[1].observables[3].expression, _ = ObservableExpression.deserialize('3 * {}'.format(models[1].species[1].id), objs)

        objs = models[2].gen_serialized_val_obj_map()
        models[2].observables[0].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[2].species[0].id), objs)
        models[2].observables[1].expression, _ = ObservableExpression.deserialize('2 * {}'.format(models[2].species[0].id), objs)
        models[2].observables[2].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[2].species[1].id), objs)
        models[2].observables[3].expression, _ = ObservableExpression.deserialize('2 * {}'.format(models[2].species[1].id), objs)
        models[2].observables[4].expression, _ = ObservableExpression.deserialize('3 * {}'.format(models[2].species[0].id), objs)
        models[2].observables[5].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[2].species[3].id), objs)
        models[2].observables[6].expression, _ = ObservableExpression.deserialize('3 * {}'.format(models[2].species[3].id), objs)

        objs = models[3].gen_serialized_val_obj_map()
        models[3].observables[0].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[3].species[0].id), objs)
        models[3].observables[1].expression, _ = ObservableExpression.deserialize('2 * {}'.format(models[3].species[0].id), objs)
        models[3].observables[2].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[3].species[1].id), objs)
        models[3].observables[3].expression, _ = ObservableExpression.deserialize('2 * {}'.format(models[3].species[1].id), objs)
        models[3].observables[4].expression, _ = ObservableExpression.deserialize('3 * {}'.format(models[3].species[0].id), objs)
        models[3].observables[5].expression, _ = ObservableExpression.deserialize('1 * {}'.format(models[3].species[3].id), objs)
        models[3].observables[6].expression, _ = ObservableExpression.deserialize('3 * {}'.format(models[3].species[3].id), objs)

        # functions
        models[0].functions.create(id='sc_0', units=unit_registry.parse_units('dimensionless'))
        models[0].functions.create(id='sc_1', units=unit_registry.parse_units('dimensionless'))

        models[1].functions.create(id='sc_0', units=unit_registry.parse_units('dimensionless'))
        models[1].functions.create(id='sc_2', units=unit_registry.parse_units('dimensionless'))

        models[2].functions.create(id='sc_0', units=unit_registry.parse_units('dimensionless'))
        models[2].functions.create(id='sc_1', units=unit_registry.parse_units('dimensionless'))
        models[2].functions.create(id='sc_2', units=unit_registry.parse_units('dimensionless'))

        models[3].functions.create(id='sc_0', units=unit_registry.parse_units('dimensionless'))
        models[3].functions.create(id='sc_1', units=unit_registry.parse_units('dimensionless'))
        models[3].functions.create(id='sc_2', units=unit_registry.parse_units('dimensionless'))

        models[0].functions[0].expression, _ = FunctionExpression.deserialize('1.', {})
        models[0].functions[1].expression, _ = FunctionExpression.deserialize('2.', {})

        models[1].functions[0].expression, _ = FunctionExpression.deserialize('1.', {})
        models[1].functions[1].expression, _ = FunctionExpression.deserialize('3.', {})

        models[2].functions[0].expression, _ = FunctionExpression.deserialize('1.', {})
        models[2].functions[1].expression, _ = FunctionExpression.deserialize('2.', {})
        models[2].functions[2].expression, _ = FunctionExpression.deserialize('3.', {})

        models[3].functions[0].expression, _ = FunctionExpression.deserialize('1.', {})
        models[3].functions[1].expression, _ = FunctionExpression.deserialize('2.', {})
        models[3].functions[2].expression, _ = FunctionExpression.deserialize('3.', {})

        # dFBA objective reactions
        models[0].dfba_obj_reactions.create(id='dfba_rxn_0', submodel=models[0].submodels[1])
        models[0].dfba_obj_reactions.create(id='dfba_rxn_1', submodel=models[0].submodels[1])

        models[1].dfba_obj_reactions.create(id='dfba_rxn_2', submodel=models[1].submodels[1])

        models[2].dfba_obj_reactions.create(id='dfba_rxn_0', submodel=models[2].submodels[1])
        models[2].dfba_obj_reactions.create(id='dfba_rxn_1', submodel=models[2].submodels[1])
        models[2].dfba_obj_reactions.create(id='dfba_rxn_2', submodel=models[2].submodels[3])

        models[3].dfba_obj_reactions.create(id='dfba_rxn_0', submodel=models[3].submodels[1])
        models[3].dfba_obj_reactions.create(id='dfba_rxn_1', submodel=models[3].submodels[1])
        models[3].dfba_obj_reactions.create(id='dfba_rxn_2', submodel=models[3].submodels[3])

        # dFBA objective species
        models[0].dfba_obj_species.create(dfba_obj_reaction=models[0].dfba_obj_reactions[0],
                                          species=models[0].species[0], value=1.)
        models[0].dfba_obj_species.create(dfba_obj_reaction=models[0].dfba_obj_reactions[0],
                                          species=models[0].species[1], value=2.)
        models[0].dfba_obj_species.create(dfba_obj_reaction=models[0].dfba_obj_reactions[1],
                                          species=models[0].species[1], value=3.)
        models[0].dfba_obj_species.create(dfba_obj_reaction=models[0].dfba_obj_reactions[1],
                                          species=models[0].species[2], value=4.)

        models[1].dfba_obj_species.create(dfba_obj_reaction=models[1].dfba_obj_reactions[0],
                                          species=models[1].species[0], value=5.)
        models[1].dfba_obj_species.create(dfba_obj_reaction=models[1].dfba_obj_reactions[0],
                                          species=models[1].species[1], value=6.)

        models[2].dfba_obj_species.create(dfba_obj_reaction=models[2].dfba_obj_reactions[0],
                                          species=models[2].species[0], value=1.)
        models[2].dfba_obj_species.create(dfba_obj_reaction=models[2].dfba_obj_reactions[0],
                                          species=models[2].species[1], value=2.)
        models[2].dfba_obj_species.create(dfba_obj_reaction=models[2].dfba_obj_reactions[1],
                                          species=models[2].species[1], value=3.)
        models[2].dfba_obj_species.create(dfba_obj_reaction=models[2].dfba_obj_reactions[1],
                                          species=models[2].species[2], value=4.)
        models[2].dfba_obj_species.create(dfba_obj_reaction=models[2].dfba_obj_reactions[2],
                                          species=models[2].species[0], value=5.)
        models[2].dfba_obj_species.create(dfba_obj_reaction=models[2].dfba_obj_reactions[2],
                                          species=models[2].species[3], value=6.)

        models[3].dfba_obj_species.create(dfba_obj_reaction=models[3].dfba_obj_reactions[0],
                                          species=models[3].species[0], value=1.)
        models[3].dfba_obj_species.create(dfba_obj_reaction=models[3].dfba_obj_reactions[0],
                                          species=models[3].species[1], value=2.)
        models[3].dfba_obj_species.create(dfba_obj_reaction=models[3].dfba_obj_reactions[1],
                                          species=models[3].species[1], value=3.)
        models[3].dfba_obj_species.create(dfba_obj_reaction=models[3].dfba_obj_reactions[1],
                                          species=models[3].species[2], value=4.)
        models[3].dfba_obj_species.create(dfba_obj_reaction=models[3].dfba_obj_reactions[2],
                                          species=models[3].species[0], value=5.)
        models[3].dfba_obj_species.create(dfba_obj_reaction=models[3].dfba_obj_reactions[2],
                                          species=models[3].species[3], value=6.)

        for model in models:
            for obj in model.dfba_obj_species:
                obj.id = obj.gen_id()

        # dFBA objectives
        models[0].dfba_objs.create(submodel=models[0].submodels[1])

        models[1].dfba_objs.create(submodel=models[1].submodels[1])

        models[2].dfba_objs.create(submodel=models[2].submodels[1])
        models[2].dfba_objs.create(submodel=models[2].submodels[3])

        models[3].dfba_objs.create(submodel=models[3].submodels[1])
        models[3].dfba_objs.create(submodel=models[3].submodels[3])

        for model in models:
            for obj in model.dfba_objs:
                obj.id = obj.gen_id()

        objs = models[0].gen_serialized_val_obj_map()
        models[0].dfba_objs[0].expression, _ = DfbaObjectiveExpression.deserialize('dfba_rxn_0 + dfba_rxn_1', objs)

        objs = models[1].gen_serialized_val_obj_map()
        models[1].dfba_objs[0].expression, _ = DfbaObjectiveExpression.deserialize('dfba_rxn_2', objs)

        objs = models[2].gen_serialized_val_obj_map()
        models[2].dfba_objs[0].expression, _ = DfbaObjectiveExpression.deserialize('dfba_rxn_0 + dfba_rxn_1', objs)
        models[2].dfba_objs[1].expression, _ = DfbaObjectiveExpression.deserialize('dfba_rxn_2', objs)

        objs = models[3].gen_serialized_val_obj_map()
        models[3].dfba_objs[0].expression, _ = DfbaObjectiveExpression.deserialize('dfba_rxn_0 + dfba_rxn_1', objs)
        models[3].dfba_objs[1].expression, _ = DfbaObjectiveExpression.deserialize('dfba_rxn_2', objs)

        # stop conditions
        models[0].stop_conditions.create(id='sc_0')
        models[0].stop_conditions.create(id='sc_1')

        models[1].stop_conditions.create(id='sc_0')
        models[1].stop_conditions.create(id='sc_2')

        models[2].stop_conditions.create(id='sc_0')
        models[2].stop_conditions.create(id='sc_1')
        models[2].stop_conditions.create(id='sc_2')

        models[3].stop_conditions.create(id='sc_0')
        models[3].stop_conditions.create(id='sc_1')
        models[3].stop_conditions.create(id='sc_2')

        models[0].stop_conditions[0].expression, _ = StopConditionExpression.deserialize('1. > 0.', {})
        models[0].stop_conditions[1].expression, _ = StopConditionExpression.deserialize('2. > 0.', {})

        models[1].stop_conditions[0].expression, _ = StopConditionExpression.deserialize('1. > 0.', {})
        models[1].stop_conditions[1].expression, _ = StopConditionExpression.deserialize('3. > 0.', {})

        models[2].stop_conditions[0].expression, _ = StopConditionExpression.deserialize('1. > 0.', {})
        models[2].stop_conditions[1].expression, _ = StopConditionExpression.deserialize('2. > 0.', {})
        models[2].stop_conditions[2].expression, _ = StopConditionExpression.deserialize('3. > 0.', {})

        models[3].stop_conditions[0].expression, _ = StopConditionExpression.deserialize('1. > 0.', {})
        models[3].stop_conditions[1].expression, _ = StopConditionExpression.deserialize('2. > 0.', {})
        models[3].stop_conditions[2].expression, _ = StopConditionExpression.deserialize('3. > 0.', {})

        # identifiers
        models[0].identifiers.create(namespace='db_0', id='id_0', submodels=models[0].submodels)
        models[0].identifiers.create(namespace='db_0', id='id_1', submodels=models[0].submodels)
        models[0].identifiers.create(namespace='db_1', id='id_0', submodels=models[0].submodels)
        models[0].identifiers.create(namespace='db_1', id='id_1', submodels=models[0].submodels)

        models[1].identifiers.create(namespace='db_0', id='id_0', submodels=models[1].submodels)
        models[1].identifiers.create(namespace='db_0', id='id_2', submodels=models[1].submodels)
        models[1].identifiers.create(namespace='db_2', id='id_0', submodels=models[1].submodels)
        models[1].identifiers.create(namespace='db_2', id='id_2', submodels=models[1].submodels)

        models[2].identifiers.create(namespace='db_0', id='id_0', submodels=models[2].submodels)
        models[2].identifiers.create(namespace='db_0', id='id_1', submodels=models[2].submodels[0:2])
        models[2].identifiers.create(namespace='db_1', id='id_0', submodels=models[2].submodels[0:2])
        models[2].identifiers.create(namespace='db_1', id='id_1', submodels=models[2].submodels[0:2])
        models[2].identifiers.create(namespace='db_0', id='id_2', submodels=models[2].submodels[2:4])
        models[2].identifiers.create(namespace='db_2', id='id_0', submodels=models[2].submodels[2:4])
        models[2].identifiers.create(namespace='db_2', id='id_2', submodels=models[2].submodels[2:4])

        models[3].identifiers.create(namespace='db_0', id='id_0', submodels=models[3].submodels)
        models[3].identifiers.create(namespace='db_0', id='id_1', submodels=models[3].submodels[0:2])
        models[3].identifiers.create(namespace='db_1', id='id_0', submodels=models[3].submodels[0:2])
        models[3].identifiers.create(namespace='db_1', id='id_1', submodels=models[3].submodels[0:2])
        models[3].identifiers.create(namespace='db_0', id='id_2', submodels=models[3].submodels[2:4])
        models[3].identifiers.create(namespace='db_2', id='id_0', submodels=models[3].submodels[2:4])
        models[3].identifiers.create(namespace='db_2', id='id_2', submodels=models[3].submodels[2:4])

        # evidence
        models[0].evidences.create(id='ev_0', submodels=models[0].submodels[0:1])
        models[0].evidences.create(id='ev_1', submodels=models[0].submodels[0:1], compartments=models[0].compartments)
        models[0].evidences.create(id='ev_2', compartments=models[0].compartments[0:1])
        models[0].evidences.create(id='ev_3', compartments=models[0].compartments[1:2])

        models[1].evidences.create(id='ev_0', submodels=models[1].submodels[0:1])
        models[1].evidences.create(id='ev_1', submodels=models[1].submodels[1:2], compartments=models[1].compartments)
        models[1].evidences.create(id='ev_4', compartments=models[1].compartments[0:1])
        models[1].evidences.create(id='ev_5', compartments=models[1].compartments[1:2])

        models[2].evidences.create(id='ev_0', submodels=models[2].submodels[0:1] + models[2].submodels[2:3])
        models[2].evidences.create(id='ev_1', submodels=models[2].submodels[0:1] + models[2].submodels[3:4],
                                   compartments=models[2].compartments)
        models[2].evidences.create(id='ev_2', compartments=models[2].compartments[0:1])
        models[2].evidences.create(id='ev_3', compartments=models[2].compartments[1:2])
        models[2].evidences.create(id='ev_4', compartments=models[2].compartments[0:1])
        models[2].evidences.create(id='ev_5', compartments=models[2].compartments[2:3])

        models[3].evidences.create(id='ev_0', submodels=models[3].submodels[0:1] + models[3].submodels[2:3])
        models[3].evidences.create(id='ev_1', submodels=models[3].submodels[0:1] + models[3].submodels[3:4],
                                   compartments=models[3].compartments)
        models[3].evidences.create(id='ev_2', compartments=models[3].compartments[0:1])
        models[3].evidences.create(id='ev_3', compartments=models[3].compartments[1:2])
        models[3].evidences.create(id='ev_4', compartments=models[3].compartments[0:1])
        models[3].evidences.create(id='ev_5', compartments=models[3].compartments[2:3])

        # references
        models[0].references.create(id='ref_0', taxon=models[0].taxon, submodels=models[0].submodels)
        models[0].references.create(id='ref_1', env=models[0].env, submodels=models[0].submodels)

        models[1].references.create(id='ref_0', taxon=models[1].taxon, submodels=models[1].submodels)
        models[1].references.create(id='ref_2', env=models[1].env, submodels=models[1].submodels)

        models[2].references.create(id='ref_0', taxon=models[2].taxon, submodels=models[2].submodels)
        models[2].references.create(id='ref_1', env=models[2].env, submodels=models[2].submodels[0:2])
        models[2].references.create(id='ref_2', env=models[2].env, submodels=models[2].submodels[2:4])

        models[3].references.create(id='ref_0', taxon=models[3].taxon, submodels=models[3].submodels)
        models[3].references.create(id='ref_1', env=models[3].env, submodels=models[3].submodels[0:2])
        models[3].references.create(id='ref_2', env=models[3].env, submodels=models[3].submodels[2:4])

        for model in models:
            for ref in model.references:
                ref.type = onto['WC:article']

        # generate ids
        for model in models:
            wc_lang.util.gen_ids(model)

        return models

    def test_gen_serialized_val_obj_map(self):
        model_1, _, _, _ = self.gen_models()
        map = model_1.gen_serialized_val_obj_map()

        self.assertEqual(map[Model], {'model_1': model_1})
        self.assertEqual(map[Taxon], {'taxon': model_1.taxon})
        self.assertEqual(map[Environment], {'env': model_1.env})
        self.assertEqual(map[Submodel], {
            'submodel_11': model_1.submodels.get_one(id='submodel_11'),
            'submodel_12': model_1.submodels.get_one(id='submodel_12'),
        })
        self.assertEqual(map[Reference], {
            'ref_0': model_1.references.get_one(id='ref_0'),
            'ref_1': model_1.references.get_one(id='ref_1'),
        })
        self.assertEqual(map[Identifier], {
            'db_0: id_0': model_1.identifiers.get_one(namespace='db_0', id='id_0'),
            'db_0: id_1': model_1.identifiers.get_one(namespace='db_0', id='id_1'),
            'db_1: id_0': model_1.identifiers.get_one(namespace='db_1', id='id_0'),
            'db_1: id_1': model_1.identifiers.get_one(namespace='db_1', id='id_1'),
        })

    def test_gen_merge_map(self):
        model_1, model_2, _, _ = self.gen_models()

        other_objs_in_self, other_objs_not_in_self = model_1.gen_merge_map(model_2)

        self.assertEqual(other_objs_in_self[model_2], model_1)
        self.assertEqual(other_objs_in_self[model_2.taxon], model_1.taxon)
        self.assertEqual(other_objs_in_self[model_2.env], model_1.env)
        self.assertEqual(other_objs_in_self[model_2.compartments.get_one(id='c_0')],
                         model_1.compartments.get_one(id='c_0'))
        self.assertEqual(other_objs_in_self[model_2.references.get_one(id='ref_0')],
                         model_1.references.get_one(id='ref_0'))
        self.assertEqual(other_objs_in_self[model_2.identifiers.get_one(namespace='db_0', id='id_0')],
                         model_1.identifiers.get_one(namespace='db_0', id='id_0'))

        self.assertIn(model_2.compartments[1], other_objs_not_in_self)
        self.assertIn(model_2.submodels[0], other_objs_not_in_self)
        self.assertIn(model_2.submodels[1], other_objs_not_in_self)
        self.assertIn(model_2.references.get_one(id='ref_2'), other_objs_not_in_self)
        self.assertIn(model_2.identifiers.get_one(namespace='db_0', id='id_2'), other_objs_not_in_self)
        self.assertIn(model_2.identifiers.get_one(namespace='db_2', id='id_0'), other_objs_not_in_self)
        self.assertIn(model_2.identifiers.get_one(namespace='db_2', id='id_2'), other_objs_not_in_self)

    def test_gen_merge_map_different_taxon_env_ids(self):
        model_1, model_2, _, _ = self.gen_models()
        model_2.taxon.id = 'taxon_2'
        model_2.env.id = 'env_2'
        other_objs_in_self, other_objs_not_in_self = model_1.gen_merge_map(model_2)

        self.assertNotIn(model_2.taxon, other_objs_in_self)
        self.assertIn(model_2.taxon, other_objs_not_in_self)

    def test_gen_merge_map_self(self):
        model_1, _, _, _ = self.gen_models()

        model_2 = model_1.copy()
        other_objs_in_self, other_objs_not_in_self = model_1.gen_merge_map(model_2)

        model_1.normalize()
        model_2.normalize()

        self.assertEqual(set(other_objs_in_self.keys()), set(model_2.get_related()))
        self.assertEqual(set(other_objs_in_self.values()), set(model_1.get_related()))
        self.assertEqual(set(model_2.get_related()), set(other_objs_in_self.keys()))
        self.assertEqual(set(model_1.get_related()), set(other_objs_in_self.values()))

        self.assertEqual(other_objs_in_self,
                         {o: s for s, o in zip(model_1.get_related(), model_2.get_related())})

        self.assertEqual(other_objs_not_in_self, [])

    def test_gen_merge_map_no_taxon_env(self):
        model_1, model_2, _, _ = self.gen_models()
        taxon_2 = model_2.taxon
        env_2 = model_2.env
        model_2.taxon.identifiers = []
        model_2.taxon.evidence = []
        model_2.taxon.references = []
        model_2.taxon = None
        model_2.env.identifiers = []
        model_2.env.evidence = []
        model_2.env.references = []
        model_2.env = None
        other_objs_in_self, other_objs_not_in_self = model_1.gen_merge_map(model_2)
        self.assertNotIn(taxon_2, other_objs_in_self)
        self.assertNotIn(env_2, other_objs_in_self)
        self.assertNotIn(taxon_2, other_objs_not_in_self)
        self.assertNotIn(env_2, other_objs_not_in_self)

        model_1, model_2, _, _ = self.gen_models()
        taxon_2 = model_2.taxon
        env_2 = model_2.env
        model_1.taxon.identifiers = []
        model_1.taxon.evidence = []
        model_1.taxon.references = []
        model_1.taxon = None
        model_1.env.identifiers = []
        model_1.env.evidence = []
        model_1.env.references = []
        model_1.env = None
        other_objs_in_self, other_objs_not_in_self = model_1.gen_merge_map(model_2)
        self.assertNotIn(taxon_2, other_objs_in_self)
        self.assertNotIn(env_2, other_objs_in_self)
        self.assertIn(taxon_2, other_objs_not_in_self)
        self.assertIn(env_2, other_objs_not_in_self)

    def test_gen_merge_map_repeated_id(self):
        model_1, model_2, _, _ = self.gen_models()
        model_1.references.get_one(id='ref_1').id = 'ref_0'
        with self.assertRaisesRegex(ValueError, 'Serialized value "ref_0" is not unique for Reference'):
            model_1.gen_merge_map(model_2)

        model_1, model_2, _, _ = self.gen_models()
        model_2.references.get_one(id='ref_2').id = 'ref_0'
        with self.assertRaisesRegex(ValueError, 'Serialized value "ref_0" is not unique for Reference'):
            model_1.gen_merge_map(model_2)

    def test_merge(self):
        model_1, model_2, model_3, _ = self.gen_models()
        model_1.merge(model_2)
        self.assertTrue(model_1.is_equal(model_3))
        self.assertNotIn(model_2, model_1.get_related())

        model_1, model_2, _, model_3 = self.gen_models()
        model_2.merge(model_1)
        self.assertTrue(model_2.is_equal(model_3))
        self.assertNotIn(model_2, model_1.get_related())

    def test_merge_taxon(self):
        model_1, model_2, _, _ = self.gen_models()

        merged_model = model_1.copy()
        model_2_copy = model_2.copy()
        merged_model.merge(model_2_copy)

        self.assertEqual(merged_model.taxon.id, model_1.taxon.id)
        self.assertEqual(merged_model.taxon.name, model_1.taxon.name)
        self.assertEqual(merged_model.taxon.rank, model_1.taxon.rank)
        self.assertEqual(merged_model.taxon.comments, model_1.taxon.comments)
        self.assertEqual(set(id.serialize() for id in merged_model.taxon.identifiers),
                         set(id.serialize() for id in model_1.taxon.identifiers) |
                         set(id.serialize() for id in model_2.taxon.identifiers))
        self.assertEqual(set(ref.id for ref in merged_model.taxon.references),
                         set(ref.id for ref in model_1.taxon.references) |
                         set(ref.id for ref in model_2.taxon.references))

        self.assertNotIn(model_2_copy, merged_model.get_related())

    def test_merge_no_secondary_taxon(self):
        model_1, model_2, _, _ = self.gen_models()
        model_2.taxon.identifiers = []
        model_2.taxon.references = []
        model_2.taxon = None

        merged_model = model_1.copy()
        model_2_copy = model_2.copy()
        merged_model.merge(model_2_copy)

        self.assertEqual(merged_model.taxon.id, model_1.taxon.id)
        self.assertEqual(merged_model.taxon.name, model_1.taxon.name)
        self.assertEqual(merged_model.taxon.rank, model_1.taxon.rank)
        self.assertEqual(merged_model.taxon.comments, model_1.taxon.comments)
        self.assertEqual(set(id.serialize() for id in merged_model.taxon.identifiers),
                         set(id.serialize() for id in model_1.taxon.identifiers))
        self.assertEqual(set(ref.id for ref in merged_model.taxon.references),
                         set(ref.id for ref in model_1.taxon.references))

        self.assertNotIn(model_2_copy, merged_model.get_related())

    def test_merge_no_primary_taxon(self):
        model_1, model_2, _, _ = self.gen_models()
        model_1.taxon.identifiers = []
        model_1.taxon.references = []
        model_1.taxon = None

        merged_model = model_1.copy()
        model_2_copy = model_2.copy()
        merged_model.merge(model_2_copy)

        self.assertEqual(merged_model.taxon.id, model_2.taxon.id)
        self.assertEqual(merged_model.taxon.name, model_2.taxon.name)
        self.assertEqual(merged_model.taxon.rank, model_2.taxon.rank)
        self.assertEqual(merged_model.taxon.comments, model_2.taxon.comments)
        self.assertEqual(set(id.serialize() for id in merged_model.taxon.identifiers),
                         set(id.serialize() for id in model_2.taxon.identifiers))
        self.assertEqual(set(ref.id for ref in merged_model.taxon.references),
                         set(ref.id for ref in model_2.taxon.references))

        self.assertNotIn(model_2_copy, merged_model.get_related())

    def test_merge_no_taxa(self):
        model_1, model_2, _, _ = self.gen_models()
        model_1.taxon.identifiers = []
        model_1.taxon.references = []
        model_1.taxon = None
        model_2.taxon.identifiers = []
        model_2.taxon.references = []
        model_2.taxon = None

        merged_model = model_1.copy()
        model_2_copy = model_2.copy()
        merged_model.merge(model_2_copy)
        self.assertEqual(merged_model.taxon, None)

        self.assertNotIn(model_2_copy, merged_model.get_related())

    def test_merge_env(self):
        model_1, model_2, _, _ = self.gen_models()

        merged_model = model_1.copy()
        model_2_copy = model_2.copy()
        merged_model.merge(model_2_copy)

        self.assertEqual(merged_model.env.id, model_1.env.id)
        self.assertEqual(merged_model.env.name, model_1.env.name)
        self.assertEqual(merged_model.env.comments, model_1.env.comments + '\n\n' + model_2.env.comments)
        self.assertEqual(set(id.serialize() for id in merged_model.env.identifiers),
                         set(id.serialize() for id in model_1.env.identifiers) |
                         set(id.serialize() for id in model_2.env.identifiers))
        self.assertEqual(set(ref.id for ref in merged_model.env.references),
                         set(ref.id for ref in model_1.env.references) |
                         set(ref.id for ref in model_2.env.references))

        self.assertNotIn(model_2_copy, merged_model.get_related())

    def test_merge_no_secondary_env(self):
        model_1, model_2, _, _ = self.gen_models()
        model_2.env.identifiers = []
        model_2.env.references = []
        model_2.env = None

        merged_model = model_1.copy()
        model_2_copy = model_2.copy()
        merged_model.merge(model_2_copy)

        self.assertEqual(merged_model.env.id, model_1.env.id)
        self.assertEqual(merged_model.env.name, model_1.env.name)
        self.assertEqual(merged_model.env.comments, model_1.env.comments)
        self.assertEqual(set(id.serialize() for id in merged_model.env.identifiers),
                         set(id.serialize() for id in model_1.env.identifiers))
        self.assertEqual(set(ref.id for ref in merged_model.env.references),
                         set(ref.id for ref in model_1.env.references))

        self.assertNotIn(model_2_copy, merged_model.get_related())

    def test_merge_no_primary_env(self):
        model_1, model_2, _, _ = self.gen_models()
        model_1.env.identifiers = []
        model_1.env.references = []
        model_1.env = None

        merged_model = model_1.copy()
        model_2_copy = model_2.copy()
        merged_model.merge(model_2_copy)

        self.assertEqual(merged_model.env.id, model_2.env.id)
        self.assertEqual(merged_model.env.name, model_2.env.name)
        self.assertEqual(merged_model.env.comments, model_2.env.comments)
        self.assertEqual(set(id.serialize() for id in merged_model.env.identifiers),
                         set(id.serialize() for id in model_2.env.identifiers))
        self.assertEqual(set(ref.id for ref in merged_model.env.references),
                         set(ref.id for ref in model_2.env.references))

        self.assertNotIn(model_2_copy, merged_model.get_related())

    def test_merge_no_env(self):
        model_1, model_2, _, _ = self.gen_models()
        model_1.env.identifiers = []
        model_1.env.references = []
        model_1.env = None
        model_2.env.identifiers = []
        model_2.env.references = []
        model_2.env = None

        merged_model = model_1.copy()
        model_2_copy = model_2.copy()
        merged_model.merge(model_2_copy)
        self.assertEqual(merged_model.env, None)

        self.assertNotIn(model_2_copy, merged_model.get_related())

    def test_merge_submodels(self):
        model_1, model_2, _, _ = self.gen_models()
        model_3 = model_1.copy()
        model_3.merge(model_2.copy())

        self.assertEqual(len(model_3.submodels), 4)

        s_11 = model_1.submodels.get_one(id='submodel_11')
        s_31 = model_3.submodels.get_one(id='submodel_11')
        self.assertEqual(set(obj.serialize() for obj in s_11.identifiers),
                         set(obj.serialize() for obj in s_31.identifiers))
        self.assertEqual(set(obj.serialize() for obj in s_11.evidence),
                         set(obj.serialize() for obj in s_31.evidence))
        self.assertEqual(set(obj.serialize() for obj in s_11.references),
                         set(obj.serialize() for obj in s_31.references))

        s_22 = model_2.submodels.get_one(id='submodel_21')
        s_32 = model_3.submodels.get_one(id='submodel_21')
        self.assertEqual(set(obj.serialize() for obj in s_22.identifiers),
                         set(obj.serialize() for obj in s_32.identifiers))
        self.assertEqual(set(obj.serialize() for obj in s_22.evidence),
                         set(obj.serialize() for obj in s_32.evidence))
        self.assertEqual(set(obj.serialize() for obj in s_22.references),
                         set(obj.serialize() for obj in s_32.references))

        self.assertEqual(s_31.identifiers.get_one(namespace='db_0', id='id_0'),
                         s_32.identifiers.get_one(namespace='db_0', id='id_0'))
        self.assertEqual(s_31.evidence.get_one(id='ref_0'),
                         s_32.evidence.get_one(id='ref_0'))
        self.assertEqual(s_31.references.get_one(id='ref_0'),
                         s_32.references.get_one(id='ref_0'))

    def test_merge_submodels_same_id(self):
        model_1, model_2, _, _ = self.gen_models()
        model_2.submodels[0].id = 'submodel_11'
        with self.assertRaisesRegex(ValueError, 'cannot be joined'):
            model_1.merge(model_2)

    def test_merge_one_to_one(self):
        model_1 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_2 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_3 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_1.compartments.create(id='c_0', geometry=None)
        model_2.compartments.create(id='c_1', geometry=None)
        model_3.compartments.create(id='c_0', geometry=None)
        model_3.compartments.create(id='c_1', geometry=None)
        model_1.parameters.create(id='p_0', density_compartment=model_1.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_2.parameters.create(id='p_1', density_compartment=model_2.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_3.parameters.create(id='p_0', density_compartment=model_3.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_3.parameters.create(id='p_1', density_compartment=model_3.compartments[1], units=unit_registry.parse_units('g l^-1'))
        merged_model = model_1.copy()
        merged_model.merge(model_2.copy())
        self.assertTrue(merged_model.is_equal(model_3))
        merged_model = model_2.copy()
        merged_model.merge(model_1.copy())
        self.assertTrue(merged_model.is_equal(model_3))

        model_1 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_2 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_3 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_1.compartments.create(id='c_0', geometry=None)
        model_2.compartments.create(id='c_1', geometry=None)
        model_3.compartments.create(id='c_0', geometry=None)
        model_3.compartments.create(id='c_1', geometry=None)
        model_1.parameters.create(id='p_0', density_compartment=model_1.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_3.parameters.create(id='p_0', density_compartment=model_3.compartments[0], units=unit_registry.parse_units('g l^-1'))
        merged_model = model_1.copy()
        merged_model.merge(model_2.copy())
        self.assertTrue(merged_model.is_equal(model_3))
        merged_model = model_2.copy()
        merged_model.merge(model_1.copy())
        self.assertTrue(merged_model.is_equal(model_3))

        model_1 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_2 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_3 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_1.compartments.create(id='c_0', geometry=None)
        model_2.compartments.create(id='c_1', geometry=None)
        model_3.compartments.create(id='c_0', geometry=None)
        model_3.compartments.create(id='c_1', geometry=None)
        model_2.parameters.create(id='p_1', density_compartment=model_2.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_3.parameters.create(id='p_1', density_compartment=model_3.compartments[1], units=unit_registry.parse_units('g l^-1'))
        merged_model = model_1.copy()
        merged_model.merge(model_2.copy())
        self.assertTrue(merged_model.is_equal(model_3))
        merged_model = model_2.copy()
        merged_model.merge(model_1.copy())
        self.assertTrue(merged_model.is_equal(model_3))

        model_1 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_2 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_3 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_1.compartments.create(id='c_0', geometry=None)
        model_2.compartments.create(id='c_1', geometry=None)
        model_3.compartments.create(id='c_0', geometry=None)
        model_3.compartments.create(id='c_1', geometry=None)
        model_1.parameters.create(id='p_1', density_compartment=model_1.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_1.parameters.create(id='p_2', units=unit_registry.parse_units('g l^-1'))
        model_2.parameters.create(id='p_1', units=unit_registry.parse_units('g l^-1'))
        model_2.parameters.create(id='p_3', units=unit_registry.parse_units('g l^-1'))
        model_3.parameters.create(id='p_1', density_compartment=model_3.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_3.parameters.create(id='p_2', units=unit_registry.parse_units('g l^-1'))
        model_3.parameters.create(id='p_3', units=unit_registry.parse_units('g l^-1'))
        merged_model = model_1.copy()
        merged_model.merge(model_2.copy())
        self.assertTrue(merged_model.is_equal(model_3))
        merged_model = model_2.copy()
        merged_model.merge(model_1.copy())
        self.assertTrue(merged_model.is_equal(model_3))

        model_1 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_2 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_3 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_1.compartments.create(id='c_0', geometry=None)
        model_2.compartments.create(id='c_0', geometry=None)
        model_3.compartments.create(id='c_0', geometry=None)
        model_1.parameters.create(id='p_1', density_compartment=model_1.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_2.parameters.create(id='p_1', units=unit_registry.parse_units('g l^-1'))
        model_3.parameters.create(id='p_1', density_compartment=model_3.compartments[0], units=unit_registry.parse_units('g l^-1'))
        merged_model = model_1.copy()
        merged_model.merge(model_2.copy())
        self.assertTrue(merged_model.is_equal(model_3))
        merged_model = model_2.copy()
        merged_model.merge(model_1.copy())
        self.assertTrue(merged_model.is_equal(model_3))

    def test_merge_one_to_one_errors(self):
        model_1 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_2 = Model(id='model_2', version='0.0.1', wc_lang_version='0.0.1')
        model_1.compartments.create(id='c_0', geometry=None)
        model_2.compartments.create(id='c_0', geometry=None)
        model_1.parameters.create(id='p_0', density_compartment=model_1.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_2.parameters.create(id='p_1', density_compartment=model_2.compartments[0], units=unit_registry.parse_units('g l^-1'))
        with self.assertRaisesRegex(ValueError, 'Cannot join'):
            model_1.copy().merge(model_2.copy())
        with self.assertRaisesRegex(ValueError, 'Cannot join'):
            model_2.copy().merge(model_1.copy())

        model_1 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_2 = Model(id='model_2', version='0.0.1', wc_lang_version='0.0.1')
        model_1.compartments.create(id='c_0', geometry=None)
        model_2.compartments.create(id='c_1', geometry=None)
        model_1.parameters.create(id='p_0', density_compartment=model_1.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_2.parameters.create(id='p_0', density_compartment=model_2.compartments[0], units=unit_registry.parse_units('g l^-1'))
        with self.assertRaisesRegex(ValueError, 'Cannot join'):
            model_1.copy().merge(model_2.copy())
        with self.assertRaisesRegex(ValueError, 'Cannot join'):
            model_2.copy().merge(model_1.copy())

        model_1 = Model(id='model_1', version='0.0.1', wc_lang_version='0.0.1')
        model_2 = Model(id='model_2', version='0.0.1', wc_lang_version='0.0.1')
        model_1.compartments.create(id='c_0', geometry=None)
        model_2.compartments.create(id='c_0', geometry=None)
        model_1.parameters.create(id='p_1', density_compartment=model_1.compartments[0], units=unit_registry.parse_units('g l^-1'))
        model_1.parameters.create(id='p_2', units=unit_registry.parse_units('g l^-1'))
        model_2.parameters.create(id='p_1', units=unit_registry.parse_units('g l^-1'))
        model_2.parameters.create(id='p_2', density_compartment=model_2.compartments[0], units=unit_registry.parse_units('g l^-1'))
        with self.assertRaisesRegex(ValueError, 'Cannot join'):
            model_1.copy().merge(model_2.copy())
        with self.assertRaisesRegex(ValueError, 'Cannot join'):
            model_2.copy().merge(model_1.copy())

    def test_cut_and_merge(self):
        _, _, model, _ = self.gen_models()

        core_model, submodels = model.submodels.gen_models()
        for submodel in submodels:
            core_model.merge(submodel)

        self.assertTrue(model.is_equal(core_model))
