""" Tests of utilities.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang.core import (Model, Taxon, Environment, Submodel,
                          Compartment,
                          SpeciesType, Species, SpeciesCoefficient, DistributionInitConcentration,
                          Reaction, RateLaw, RateLawExpression, Parameter,
                          DfbaObjSpecies, DfbaObjReaction,
                          DfbaObjective, DfbaObjectiveExpression,
                          Observable, ObservableExpression,
                          Function, FunctionExpression,
                          StopCondition, StopConditionExpression,
                          Observation, ObservationSet, Evidence, Conclusion,
                          Reference, Author, Change, Identifier,
                          InitVolume, Ph, ChemicalStructure, FluxBounds, ObservationGenotype, ObservationEnv, Process,
                          )
from wc_lang import io
from wc_lang import util
from wc_onto import onto
from wc_utils.util.units import unit_registry
import obj_tables.sci.units
import os.path
import shutil
import tempfile
import unittest


class TestUtil(unittest.TestCase):
    """ Test utilities """

    def setUp(self):
        self.model = mdl = Model()

        self.comp_0 = comp_0 = mdl.compartments.create(id='comp_0', name='compartment 0')
        self.comp_1 = comp_1 = mdl.compartments.create(id='comp_1', name='compartment 1')
        self.compartments = compartments = [comp_0, comp_1]

        self.species_types = species_types = []
        self.species = species = []
        for i in range(8):
            spec_type = mdl.species_types.create(id='spec_type_{}'.format(
                i), name='species type {}'.format(i), type=onto['WC:metabolite'])
            species_types.append(spec_type)

            if i != 3:
                spec = Species(species_type=spec_type, compartment=comp_0)
            else:
                spec = Species(species_type=spec_type, compartment=comp_1)
            spec.id = spec.gen_id()
            spec.model = mdl
            species.append(spec)

            conc = DistributionInitConcentration(species=spec, mean=1)
            conc.id = conc.gen_id()
            conc.model = mdl

        self.submdl_0 = submdl_0 = mdl.submodels.create(id='submdl_0', framework=onto['WC:stochastic_simulation_algorithm'])
        self.submdl_1 = submdl_1 = mdl.submodels.create(id='submdl_1', framework=onto['WC:stochastic_simulation_algorithm'])
        self.submdl_2 = submdl_2 = mdl.submodels.create(id='submdl_2', framework=onto['WC:dynamic_flux_balance_analysis'])
        self.submodels = [submdl_0, submdl_1, submdl_2]

        self.rxn_0 = rxn_0 = submdl_0.reactions.create(id='rxn_0', model=mdl)
        rxn_0.participants.create(species=species[0], coefficient=-2)
        rxn_0.participants.create(species=species[1], coefficient=-3)
        rxn_0.participants.create(species=species[2], coefficient=1)
        expression = RateLawExpression(
            expression='k_cat_0 * {0} / (k_m_0 + {0})'.format(species[5].get_primary_attribute()),
            species=species[5:6])
        expression.parameters.create(id='k_cat_0', value=2, model=mdl)
        expression.parameters.create(id='k_m_0', value=1, model=mdl)
        rate_law_0 = rxn_0.rate_laws.create(expression=expression, model=mdl)

        self.rxn_1 = rxn_1 = submdl_1.reactions.create(id='rxn_1', model=mdl)
        rxn_1.participants.create(species=species[0], coefficient=-2)
        rxn_1.participants.create(species=species[1], coefficient=-3)
        rxn_1.participants.create(species=species[3], coefficient=2)
        expression = RateLawExpression(
            expression='k_cat_1 * {0} / (k_m_1 + {0})'.format(species[6].get_primary_attribute()),
            species=species[6:7])
        expression.parameters.create(id='k_cat_1', value=2, model=mdl)
        expression.parameters.create(id='k_m_1', value=1, model=mdl)
        rate_law_1 = rxn_1.rate_laws.create(expression=expression, model=mdl)

        self.rxn_2 = rxn_2 = submdl_2.reactions.create(id='rxn_2', model=mdl)
        rxn_2.participants.create(species=species[0], coefficient=-2)
        rxn_2.participants.create(species=species[1], coefficient=-3)
        rxn_2.participants.create(species=species[4], coefficient=1)
        expression = RateLawExpression(
            expression='k_cat_2 * {0} / (k_m_2 + {0})'.format(species[7].get_primary_attribute()),
            species=species[7:8])
        expression.parameters.create(id='k_cat_2', value=2, model=mdl)
        expression.parameters.create(id='k_m_2', value=1, model=mdl)
        rate_law_2 = rxn_2.rate_laws.create(expression=expression, model=mdl)

        self.reactions = [rxn_0, rxn_1, rxn_2]
        self.rate_laws = [rate_law_0, rate_law_1, rate_law_2]

        self.parameters = parameters = []
        self.references = references = []
        self.identifiers = identifiers = []
        for i in range(3):
            param = mdl.parameters.create(id='param_{}'.format(i))
            parameters.append(param)

            ref = param.references.create(id='ref_{}'.format(i), type=None)
            ref.model = mdl
            references.append(ref)

            x_ref = ref.identifiers.create(namespace='Y', id='x')
            identifiers.append(x_ref)

    def test_get_model_size(self):
        model = self.model
        size = util.get_model_size(model)
        self.assertEqual(3, size['submodels'])
        self.assertEqual(8, size['species_types'])
        self.assertEqual(8, size['species'])
        self.assertEqual(3, size['reactions'])
        self.assertEqual(2, size['compartments'])
        self.assertEqual(9, size['parameters'])
        self.assertEqual(3, size['references'])

    def test_get_model_summary(self):
        model = self.model
        summary = util.get_model_summary(model)
        self.assertIsInstance(summary, str)

    def test_get_models(self):
        non_inline_models = set([
            Model, Taxon, Environment,
            Submodel, Compartment, SpeciesType, Species, Observable, DistributionInitConcentration,
            DfbaObjective,
            Reaction, RateLaw, DfbaObjSpecies, DfbaObjReaction, Parameter, Function,
            StopCondition, Observation, ObservationSet, Conclusion, Reference, Author, Change,
        ])
        inline_models = set([
            SpeciesCoefficient, RateLawExpression,
            DfbaObjectiveExpression, FunctionExpression, StopConditionExpression, ObservableExpression,
            Identifier, InitVolume, Ph, ChemicalStructure, FluxBounds, ObservationGenotype, ObservationEnv, Process,
            Evidence,
        ])
        self.assertEqual(set(util.get_models()), non_inline_models | inline_models)
        self.assertEqual(set(util.get_models(inline=False)), non_inline_models)

    def test_gen_ids(self):
        model = Model()
        model.compartments.create(id='c_1')
        model.compartments.create(id='c_2')
        model.species_types.create(id='st_1')
        model.species_types.create(id='st_2')
        model.species.create(species_type=model.species_types[0], compartment=model.compartments[0])
        model.species.create(species_type=model.species_types[0], compartment=model.compartments[1])
        model.species.create(species_type=model.species_types[1], compartment=model.compartments[0])
        model.species.create(species_type=model.species_types[1], compartment=model.compartments[1])
        util.gen_ids(model)
        self.assertEqual(model.species[0].id, 'st_1[c_1]')
        self.assertEqual(model.species[1].id, 'st_1[c_2]')
        self.assertEqual(model.species[2].id, 'st_2[c_1]')
        self.assertEqual(model.species[3].id, 'st_2[c_2]')

    def test_get_obj_units(self):
        model = Model()
        units = set([model.time_units])
        self.assertEqual(set(obj_tables.sci.units.get_obj_units(model)), units)

        model.compartments.create()
        model.compartments.create()
        for c in model.compartments:
            units.add(c.mass_units)
            if c.init_volume:
                units.add(c.init_volume.units)
            if c.ph:
                units.add(c.ph.units)
        self.assertEqual(set(obj_tables.sci.units.get_obj_units(model)), units)

        model.species_types.create()
        model.species_types.create()
        self.assertEqual(set(obj_tables.sci.units.get_obj_units(model)), units)

        for c in model.compartments:
            for s in model.species_types:
                model.species.create(compartment=c, species_type=s)
        for s in model.species:
            units.add(s.units)
        self.assertEqual(set(obj_tables.sci.units.get_obj_units(model)), units)

        model.distribution_init_concentrations.create(species=model.species[0], units=unit_registry.parse_units('M'))
        model.distribution_init_concentrations.create(species=model.species[1], units=unit_registry.parse_units('molecule'))
        for o in model.distribution_init_concentrations:
            units.add(o.units)
        self.assertEqual(set(obj_tables.sci.units.get_obj_units(model)), units)

        model.parameters.create(units=unit_registry.parse_units('g'))
        model.parameters.create(units=unit_registry.parse_units('l'))
        model.parameters.create(units=unit_registry.parse_units('s'))
        for p in model.parameters:
            units.add(p.units)
        self.assertEqual(set(obj_tables.sci.units.get_obj_units(model)), units)

        model.functions.create(units=unit_registry.parse_units('g / l'))
        model.functions.create(units=unit_registry.parse_units('g / s'))
        model.functions.create(units=unit_registry.parse_units('l / s'))
        for f in model.functions:
            units.add(f.units)
        self.assertEqual(set(obj_tables.sci.units.get_obj_units(model)), units)
