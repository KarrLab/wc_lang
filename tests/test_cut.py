""" Tests of cutting

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-28
:Copyright: 2019, Karr Lab
:License: MIT
"""

import datetime
import unittest
from wc_lang.core import Model, Species, Parameter, RateLawExpression, Reference, StopConditionExpression, ChemicalStructure
from wc_utils.util.units import unit_registry


def gen_model(submodels=True, extra_species=True):
    timestamp = datetime.datetime(2018, 1, 1, 12, 0, 0)
    model = Model(id='model', version='0.0.1', wc_lang_version='0.0.2', created=timestamp, updated=timestamp)

    if submodels:
        submodel_0 = model.submodels.create(id='submodel_0')
        submodel_1 = model.submodels.create(id='submodel_1')
    else:
        submodel_0 = None
        submodel_1 = None

    model.compartments.create(id='c_0')
    model.compartments.create(id='c_1')
    structure = ChemicalStructure(charge=0)
    model.species_types.create(id='s_0', structure=structure)
    model.species_types.create(id='s_1', structure=structure)
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[0])
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[1])
    model.species.create(compartment=model.compartments[1], species_type=model.species_types[0])
    if extra_species:
        model.species.create(compartment=model.compartments[1], species_type=model.species_types[1])
    for species in model.species:
        species.id = species.gen_id()

    if submodels:
        model.reactions.create(id='rxn_00', submodel=submodel_0)
        model.reactions.create(id='rxn_01', submodel=submodel_0)
        model.reactions.create(id='rxn_10', submodel=submodel_1)
        model.reactions.create(id='rxn_11', submodel=submodel_1)
        model.reactions[0].participants.create(species=model.species[0], coefficient=-1)
        model.reactions[0].participants.create(species=model.species[1], coefficient=1)
        model.reactions[1].participants.add(model.reactions[0].participants[0])
        model.reactions[2].participants.add(model.reactions[0].participants[0])
        model.reactions[2].participants.create(species=model.species[2], coefficient=1)
        model.reactions[3].participants.add(model.reactions[0].participants[0])

    model.parameters.create(id='p_0', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_1', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_2', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_3', value=1., units=unit_registry.parse_units('g l^-1'))
    model.parameters.create(id='p_4', value=1., units=unit_registry.parse_units('g l^-1'))
    model.compartments[0].init_density = model.parameters[3]
    model.compartments[1].init_density = model.parameters[4]
    model.parameters.create(id='p_5', value=1., units=unit_registry.parse_units('s^-1'))
    model.parameters.create(id='p_6', value=1., units=unit_registry.parse_units('s^-1'))
    model.parameters.create(id='p_7', value=1., units=unit_registry.parse_units('s^-1'))
    model.parameters.create(id='p_8', value=1., units=unit_registry.parse_units('s^-1'))
    if submodels:
        rl = model.rate_laws.create(reaction=model.reactions[0])
        rl.expression, error = RateLawExpression.deserialize(f'{model.parameters[5].id}',
                                                             {Parameter: {model.parameters[5].id: model.parameters[5]}})
        rl = model.rate_laws.create(reaction=model.reactions[1])
        rl.expression, error = RateLawExpression.deserialize(f'{model.parameters[6].id}',
                                                             {Parameter: {model.parameters[6].id: model.parameters[6]}})
        rl = model.rate_laws.create(reaction=model.reactions[2])
        rl.expression, error = RateLawExpression.deserialize(f'{model.parameters[7].id}',
                                                             {Parameter: {model.parameters[7].id: model.parameters[7]}})
        rl = model.rate_laws.create(reaction=model.reactions[3])
        rl.expression, error = RateLawExpression.deserialize(f'{model.parameters[8].id}',
                                                             {Parameter: {model.parameters[8].id: model.parameters[8]}})
        for rl in model.rate_laws:
            rl.id = rl.gen_id()

    objs = {
        Species: {s.id: s for s in model.species},
        Parameter: {p.id: p for p in model.parameters},
    }
    model.stop_conditions.create(id='sc_0')
    model.stop_conditions.create(id='sc_1')
    model.stop_conditions.create(id='sc_2')
    model.stop_conditions[0].expression, error = StopConditionExpression.deserialize(
        f'{model.species[0].id} > {model.parameters[0].id}', objs)
    assert error is None, str(error)
    model.stop_conditions[1].expression, error = StopConditionExpression.deserialize(
        f'{model.species[1].id} > {model.parameters[1].id}', objs)
    assert error is None, str(error)
    model.stop_conditions[2].expression, error = StopConditionExpression.deserialize(
        f'{model.species[2].id} > {model.parameters[2].id}', objs)
    assert error is None, str(error)

    if submodels:
        model.references.create(id='ref_0', submodels=[submodel_0, submodel_1])
        model.references.create(id='ref_1', submodels=[submodel_0])
        model.references.create(id='ref_2', submodels=[submodel_1])

    return (model, submodel_0, submodel_1)


def gen_core_model(extra_species=True):
    timestamp = datetime.datetime(2018, 1, 1, 12, 0, 0)
    model = Model(id='model', version='0.0.1', wc_lang_version='0.0.2', created=timestamp, updated=timestamp)

    model.compartments.create(id='c_0')
    model.compartments.create(id='c_1')
    structure = ChemicalStructure(charge=0)
    model.species_types.create(id='s_0', structure=structure)
    model.species_types.create(id='s_1', structure=structure)
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[0])
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[1])
    model.species.create(compartment=model.compartments[1], species_type=model.species_types[0])
    if extra_species:
        model.species.create(compartment=model.compartments[1], species_type=model.species_types[1])
    for species in model.species:
        species.id = species.gen_id()

    model.parameters.create(id='p_0', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_1', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_2', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_3', value=1., units=unit_registry.parse_units('g l^-1'))
    model.parameters.create(id='p_4', value=1., units=unit_registry.parse_units('g l^-1'))
    model.compartments[0].init_density = model.parameters[3]
    model.compartments[1].init_density = model.parameters[4]

    objs = {
        Species: {s.id: s for s in model.species},
        Parameter: {p.id: p for p in model.parameters},
    }
    model.stop_conditions.create(id='sc_0')
    model.stop_conditions.create(id='sc_1')
    model.stop_conditions.create(id='sc_2')
    model.stop_conditions[0].expression, error = StopConditionExpression.deserialize(
        f'{model.species[0].id} > {model.parameters[0].id}', objs)
    assert error is None, str(error)
    model.stop_conditions[1].expression, error = StopConditionExpression.deserialize(
        f'{model.species[1].id} > {model.parameters[1].id}', objs)
    assert error is None, str(error)
    model.stop_conditions[2].expression, error = StopConditionExpression.deserialize(
        f'{model.species[2].id} > {model.parameters[2].id}', objs)
    assert error is None, str(error)

    return model


def gen_submodel_0():
    timestamp = datetime.datetime(2018, 1, 1, 12, 0, 0)
    model = Model(id='model', version='0.0.1', wc_lang_version='0.0.2', created=timestamp, updated=timestamp)

    submodel_0 = model.submodels.create(id='submodel_0')

    model.compartments.create(id='c_0')
    structure = ChemicalStructure(charge=0)
    model.species_types.create(id='s_0', structure=structure)
    model.species_types.create(id='s_1', structure=structure)
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[0])
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[1])
    for species in model.species:
        species.id = species.gen_id()

    model.reactions.create(id='rxn_00', submodel=submodel_0)
    model.reactions.create(id='rxn_01', submodel=submodel_0)
    model.reactions[0].participants.create(species=model.species[0], coefficient=-1)
    model.reactions[0].participants.create(species=model.species[1], coefficient=1)
    model.reactions[1].participants.add(model.reactions[0].participants[0])

    model.parameters.create(id='p_0', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_1', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_3', value=1., units=unit_registry.parse_units('g l^-1'))
    model.compartments[0].init_density = model.parameters[2]
    model.parameters.create(id='p_5', value=1., units=unit_registry.parse_units('s^-1'))
    model.parameters.create(id='p_6', value=1., units=unit_registry.parse_units('s^-1'))
    rl = model.rate_laws.create(reaction=model.reactions[0])
    rl.expression, error = RateLawExpression.deserialize(f'{model.parameters[3].id}',
                                                         {Parameter: {model.parameters[3].id: model.parameters[3]}})
    rl = model.rate_laws.create(reaction=model.reactions[1])
    rl.expression, error = RateLawExpression.deserialize(f'{model.parameters[4].id}',
                                                         {Parameter: {model.parameters[4].id: model.parameters[4]}})
    for rl in model.rate_laws:
        rl.id = rl.gen_id()

    objs = {
        Species: {s.id: s for s in model.species},
        Parameter: {p.id: p for p in model.parameters},
    }
    model.stop_conditions.create(id='sc_0')
    model.stop_conditions.create(id='sc_1')
    model.stop_conditions[0].expression, error = StopConditionExpression.deserialize(
        f'{model.species[0].id} > {model.parameters[0].id}', objs)
    assert error is None, str(error)
    model.stop_conditions[1].expression, error = StopConditionExpression.deserialize(
        f'{model.species[1].id} > {model.parameters[1].id}', objs)
    assert error is None, str(error)

    model.references.create(id='ref_0', submodels=[submodel_0])
    model.references.create(id='ref_1', submodels=[submodel_0])

    return model


def gen_submodel_1():
    timestamp = datetime.datetime(2018, 1, 1, 12, 0, 0)
    model = Model(id='model', version='0.0.1', wc_lang_version='0.0.2', created=timestamp, updated=timestamp)

    submodel_1 = model.submodels.create(id='submodel_1')

    model.compartments.create(id='c_0')
    model.compartments.create(id='c_1')
    structure = ChemicalStructure(charge=0)
    model.species_types.create(id='s_0', structure=structure)
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[0])
    model.species.create(compartment=model.compartments[1], species_type=model.species_types[0])
    for species in model.species:
        species.id = species.gen_id()

    model.reactions.create(id='rxn_10', submodel=submodel_1)
    model.reactions.create(id='rxn_11', submodel=submodel_1)
    model.reactions[0].participants.create(species=model.species[0], coefficient=-1)
    model.reactions[0].participants.create(species=model.species[1], coefficient=1)
    model.reactions[1].participants.add(model.reactions[0].participants[0])

    model.parameters.create(id='p_0', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_2', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_3', value=1., units=unit_registry.parse_units('g l^-1'))
    model.parameters.create(id='p_4', value=1., units=unit_registry.parse_units('g l^-1'))
    model.compartments[0].init_density = model.parameters[2]
    model.compartments[1].init_density = model.parameters[3]
    model.parameters.create(id='p_7', value=1., units=unit_registry.parse_units('s^-1'))
    model.parameters.create(id='p_8', value=1., units=unit_registry.parse_units('s^-1'))
    rl = model.rate_laws.create(reaction=model.reactions[0])
    rl.expression, error = RateLawExpression.deserialize(f'{model.parameters[4].id}',
                                                         {Parameter: {model.parameters[4].id: model.parameters[4]}})
    rl = model.rate_laws.create(reaction=model.reactions[1])
    rl.expression, error = RateLawExpression.deserialize(f'{model.parameters[5].id}',
                                                         {Parameter: {model.parameters[5].id: model.parameters[5]}})
    for rl in model.rate_laws:
        rl.id = rl.gen_id()

    objs = {
        Species: {s.id: s for s in model.species},
        Parameter: {p.id: p for p in model.parameters},
    }
    model.stop_conditions.create(id='sc_0')
    model.stop_conditions.create(id='sc_2')
    model.stop_conditions[0].expression, error = StopConditionExpression.deserialize(
        f'{model.species[0].id} > {model.parameters[0].id}', objs)
    assert error is None, str(error)
    model.stop_conditions[1].expression, error = StopConditionExpression.deserialize(
        f'{model.species[1].id} > {model.parameters[1].id}', objs)
    assert error is None, str(error)

    model.references.create(id='ref_0', submodels=[submodel_1])
    model.references.create(id='ref_2', submodels=[submodel_1])

    return model


class GetChildrenTestCase(unittest.TestCase):
    def test_Submodel_get_immediate_children(self):
        model, submodel_0, submodel_1 = gen_model()

        self.assertEqual(set(submodel_0.get_immediate_children(kind='submodel', __include_stop_conditions=False)),
                         set([model] + model.reactions[0:2] + model.references[0:2]))
        self.assertEqual(set(submodel_1.get_immediate_children(kind='submodel', __include_stop_conditions=False)),
                         set([model] + model.reactions[2:4] + [model.references[0], model.references[2]]))

        self.assertEqual(set(submodel_0.get_immediate_children(kind='submodel')),
                         set([model] + model.reactions[0:2] + model.references[0:2] + model.stop_conditions[0:2]))
        self.assertEqual(set(submodel_1.get_immediate_children(kind='submodel')),
                         set([model] + model.reactions[2:4] + [model.references[0], model.references[2],
                                                               model.stop_conditions[0], model.stop_conditions[2]]))

        self.assertEqual(set(submodel_0.get_immediate_children(kind='submodel', __type=Reference)),
                         set(model.references[0:2]))
        self.assertEqual(set(submodel_1.get_immediate_children(kind='submodel', __type=Reference)),
                         set([model.references[0], model.references[2]]))

        self.assertEqual(submodel_0.get_immediate_children(kind='submodel', id='rxn_01'),
                         model.reactions[1:2])
        self.assertEqual(submodel_1.get_immediate_children(kind='submodel', id='model'),
                         [model])

    def test_Submodel_get_children_immediate_only(self):
        model, submodel_0, submodel_1 = gen_model()

        self.assertEqual(set(submodel_0.get_children(recursive=False, kind='submodel', __include_stop_conditions=False)),
                         set([model] + model.reactions[0:2] + model.references[0:2]))
        self.assertEqual(set(submodel_1.get_children(recursive=False, kind='submodel', __include_stop_conditions=False)),
                         set([model] + model.reactions[2:4] + [model.references[0], model.references[2]]))

        self.assertEqual(set(submodel_0.get_children(recursive=False, kind='submodel')),
                         set([model] + model.reactions[0:2] + model.references[0:2] + model.stop_conditions[0:2]))
        self.assertEqual(set(submodel_1.get_children(recursive=False, kind='submodel')),
                         set([model] + model.reactions[2:4] + [model.references[0], model.references[2],
                                                               model.stop_conditions[0], model.stop_conditions[2]]))

        self.assertEqual(set(submodel_0.get_children(recursive=False, kind='submodel', __type=Reference)),
                         set(model.references[0:2]))
        self.assertEqual(set(submodel_1.get_children(recursive=False, kind='submodel', __type=Reference)),
                         set([model.references[0], model.references[2]]))

        self.assertEqual(submodel_0.get_children(recursive=False, kind='submodel', id='rxn_01'),
                         model.reactions[1:2])
        self.assertEqual(submodel_1.get_children(recursive=False, kind='submodel', id='model'),
                         [model])

    def test_Submodel_get_children(self):
        model, submodel_0, submodel_1 = gen_model()

        self.assertEqual(set(submodel_0.get_children(kind='submodel', __include_stop_conditions=False)),
                         set([model]
                             + model.reactions[0:2]
                             + model.references[0:2]
                             + model.compartments[0:1]
                             + model.species_types[0:2]
                             + [st.structure for st in model.species_types[0:2]]
                             + model.species[0:2]
                             + model.species[0].species_coefficients[0:1]
                             + model.species[1].species_coefficients[0:1]                             
                             + model.rate_laws[0:2]
                             + [model.rate_laws[0].expression, model.rate_laws[1].expression]
                             + [model.parameters[3], model.parameters[5], model.parameters[6]]
                             ))
        self.assertEqual(set(submodel_1.get_children(kind='submodel', __include_stop_conditions=False)),
                         set([model]
                             + model.reactions[2:4]
                             + [model.references[0], model.references[2]]
                             + model.compartments[0:2]
                             + model.species_types[0:1]
                             + [st.structure for st in model.species_types[0:1]]
                             + [model.species[0], model.species[2]]
                             + model.species[0].species_coefficients[0:1]
                             + model.species[2].species_coefficients[0:1]
                             + model.rate_laws[2:4]
                             + [model.rate_laws[2].expression, model.rate_laws[3].expression]
                             + model.parameters[3:5] + model.parameters[7:9]
                             ))

        self.assertEqual(set(submodel_0.get_children(kind='submodel')),
                         set([model]
                             + model.reactions[0:2]
                             + model.references[0:2]
                             + model.stop_conditions[0:2]
                             + [model.stop_conditions[0].expression, model.stop_conditions[1].expression]
                             + model.parameters[0:2]
                             + model.compartments[0:1]
                             + model.species_types[0:2]
                             + [st.structure for st in model.species_types[0:2]]
                             + model.species[0:2]
                             + model.species[0].species_coefficients[0:1]
                             + model.species[1].species_coefficients[0:1]
                             + model.rate_laws[0:2]
                             + [model.rate_laws[0].expression, model.rate_laws[1].expression]
                             + [model.parameters[3], model.parameters[5], model.parameters[6]]
                             ))
        self.assertEqual(set(submodel_1.get_children(kind='submodel')),
                         set([model]
                             + model.reactions[2:4]
                             + [model.references[0], model.references[2]]
                             + [model.stop_conditions[0], model.stop_conditions[2]]
                             + [model.stop_conditions[0].expression, model.stop_conditions[2].expression]
                             + [model.parameters[0], model.parameters[2]]
                             + model.compartments[0:2]
                             + model.species_types[0:1]
                             + [st.structure for st in model.species_types[0:1]]
                             + [model.species[0], model.species[2]]
                             + model.species[0].species_coefficients[0:1]
                             + model.species[2].species_coefficients[0:1]
                             + model.rate_laws[2:4]
                             + [model.rate_laws[2].expression, model.rate_laws[3].expression]
                             + model.parameters[3:5] + model.parameters[7:9]
                             ))

        self.assertEqual(set(submodel_0.get_children(kind='submodel', __type=Reference)),
                         set(model.references[0:2]))
        self.assertEqual(set(submodel_1.get_children(kind='submodel', __type=Reference)),
                         set([model.references[0], model.references[2]]))

        self.assertEqual(set(submodel_0.get_children(kind='submodel', __type=Parameter)),
                         set(model.parameters[0:2] + [model.parameters[3], model.parameters[5], model.parameters[6]]))
        self.assertEqual(set(submodel_1.get_children(kind='submodel', __type=Parameter)),
                         set([model.parameters[0], model.parameters[2]] + model.parameters[3:5] + model.parameters[7:9]))

        self.assertEqual(submodel_0.get_children(kind='submodel', id='rxn_01', __check_attr_defined=False),
                         model.reactions[1:2])
        self.assertEqual(submodel_1.get_children(kind='submodel', id='model', __check_attr_defined=False),
                         [model])

        self.assertEqual(submodel_0.get_children(kind='submodel', id='p_0', __check_attr_defined=False),
                         model.parameters[0:1])
        self.assertEqual(submodel_1.get_children(kind='submodel', id='p_2', __check_attr_defined=False),
                         model.parameters[2:3])


class CutTestCase(unittest.TestCase):
    def test_cut(self):
        model, _, _ = gen_model()
        submodels = model.submodels.cut(kind='submodel')

        submodel_0 = next(s for s in submodels if s.id == 'submodel_0')
        submodel_1 = next(s for s in submodels if s.id == 'submodel_1')

        self.assertTrue(submodel_0.model.is_equal(gen_submodel_0()))
        self.assertTrue(submodel_1.model.is_equal(gen_submodel_1()))

    def test_gen_models_with_1_submodel(self):
        model = Model(id='model')
        model.submodels.create(id='submodel')

        core, submodels = model.submodels.gen_models()
        self.assertTrue(model.is_equal(core))
        self.assertEqual(submodels, [])

    def test_gen_models(self):
        model, _, _ = gen_model()
        core_model, submodels = model.submodels.gen_models()
        submodel_0 = next(s for s in submodels if s.submodels[0].id == 'submodel_0')
        submodel_1 = next(s for s in submodels if s.submodels[0].id == 'submodel_1')

        self.assertTrue(submodel_0.is_equal(gen_submodel_0()))
        self.assertTrue(submodel_1.is_equal(gen_submodel_1()))
        self.assertTrue(core_model.is_equal(gen_core_model()))

        model, _, _ = gen_model(extra_species=False)
        core_model, submodels = model.submodels.gen_models()
        submodel_0 = next(s for s in submodels if s.submodels[0].id == 'submodel_0')
        submodel_1 = next(s for s in submodels if s.submodels[0].id == 'submodel_1')
        self.assertTrue(submodel_0.is_equal(gen_submodel_0()))
        self.assertTrue(submodel_1.is_equal(gen_submodel_1()))
        self.assertTrue(core_model.is_equal(gen_core_model(extra_species=False)))

        model, _, _ = gen_model(submodels=False)
        core_model, submodels = model.submodels.gen_models()
        self.assertEqual(submodels, [])
        self.assertTrue(core_model.is_equal(model))

    def test_cut_and_merge(self):
        model, _, _ = gen_model()

        core_model, submodels = model.submodels.gen_models()
        for submodel in submodels:
            core_model.merge(submodel)

        self.assertTrue(model.is_equal(core_model))
