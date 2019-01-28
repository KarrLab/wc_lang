""" Tests of cutting

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-28
:Copyright: 2019, Karr Lab
:License: MIT
"""

import unittest
from wc_lang.core import Model, Species, Parameter, Reference, StopConditionExpression
from wc_utils.util.units import unit_registry


def gen_model():
    model = Model(id='model')

    submodel_0 = model.submodels.create(id='submodel_0')
    submodel_1 = model.submodels.create(id='submodel_1')

    model.compartments.create(id='c_0')
    model.compartments.create(id='c_1')
    model.species_types.create(id='s_0')
    model.species_types.create(id='s_1')
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[0])
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[1])
    model.species.create(compartment=model.compartments[1], species_type=model.species_types[0])
    model.species.create(compartment=model.compartments[1], species_type=model.species_types[1])
    for species in model.species:
        species.id = species.gen_id()

    model.reactions.create(id='rxn_00', submodel=submodel_0)
    model.reactions.create(id='rxn_01', submodel=submodel_0)
    model.reactions.create(id='rxn_10', submodel=submodel_1)
    model.reactions.create(id='rxn_11', submodel=submodel_1)
    model.reactions[0].participants.create(species=model.species[0], coefficient=-1)
    model.reactions[0].participants.create(species=model.species[1], coefficient=1)
    model.reactions[2].participants.add(model.reactions[0].participants[0])
    model.reactions[2].participants.create(species=model.species[2], coefficient=1)

    model.parameters.create(id='p_0', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_1', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_2', value=1., units=unit_registry.parse_units('molecule'))

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

    model.references.create(id='ref_0', submodels=[submodel_0, submodel_1])
    model.references.create(id='ref_1', submodels=[submodel_0])
    model.references.create(id='ref_2', submodels=[submodel_1])

    return (model, submodel_0, submodel_1)


def gen_submodel_0():
    model = Model(id='model')

    submodel_0 = model.submodels.create(id='submodel_0')

    model.compartments.create(id='c_0')
    model.species_types.create(id='s_0')
    model.species_types.create(id='s_1')
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[0])
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[1])
    for species in model.species:
        species.id = species.gen_id()

    model.reactions.create(id='rxn_00', submodel=submodel_0)
    model.reactions.create(id='rxn_01', submodel=submodel_0)
    model.reactions[0].participants.create(species=model.species[0], coefficient=-1)
    model.reactions[0].participants.create(species=model.species[1], coefficient=1)

    model.parameters.create(id='p_0', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_1', value=1., units=unit_registry.parse_units('molecule'))

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

    return submodel_0


def gen_submodel_1():
    model = Model(id='model')

    submodel_1 = model.submodels.create(id='submodel_1')

    model.compartments.create(id='c_0')
    model.compartments.create(id='c_1')
    model.species_types.create(id='s_0')
    model.species.create(compartment=model.compartments[0], species_type=model.species_types[0])
    model.species.create(compartment=model.compartments[1], species_type=model.species_types[0])
    for species in model.species:
        species.id = species.gen_id()

    model.reactions.create(id='rxn_10', submodel=submodel_1)
    model.reactions.create(id='rxn_11', submodel=submodel_1)
    model.reactions[0].participants.create(species=model.species[0], coefficient=-1)
    model.reactions[0].participants.create(species=model.species[1], coefficient=1)

    model.parameters.create(id='p_0', value=1., units=unit_registry.parse_units('molecule'))
    model.parameters.create(id='p_2', value=1., units=unit_registry.parse_units('molecule'))

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

    return submodel_1


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
                             + model.species[0:2]
                             + model.species[0].species_coefficients[0:1]
                             + model.species[1].species_coefficients[0:1]
                             ))
        self.assertEqual(set(submodel_1.get_children(kind='submodel', __include_stop_conditions=False)),
                         set([model]
                             + model.reactions[2:4]
                             + [model.references[0], model.references[2]]
                             + model.compartments[0:2]
                             + model.species_types[0:1]
                             + [model.species[0], model.species[2]]
                             + model.species[0].species_coefficients[0:1]
                             + model.species[2].species_coefficients[0:1]))

        self.assertEqual(set(submodel_0.get_children(kind='submodel')),
                         set([model]
                             + model.reactions[0:2]
                             + model.references[0:2]
                             + model.stop_conditions[0:2]
                             + [model.stop_conditions[0].expression, model.stop_conditions[1].expression]
                             + model.parameters[0:2]
                             + model.compartments[0:1]
                             + model.species_types[0:2]
                             + model.species[0:2]
                             + model.species[0].species_coefficients[0:1]
                             + model.species[1].species_coefficients[0:1]))
        self.assertEqual(set(submodel_1.get_children(kind='submodel')),
                         set([model]
                             + model.reactions[2:4]
                             + [model.references[0], model.references[2]]
                             + [model.stop_conditions[0], model.stop_conditions[2]]
                             + [model.stop_conditions[0].expression, model.stop_conditions[2].expression]
                             + [model.parameters[0], model.parameters[2]]
                             + model.compartments[0:2]
                             + model.species_types[0:1]
                             + [model.species[0], model.species[2]]
                             + model.species[0].species_coefficients[0:1]
                             + model.species[2].species_coefficients[0:1]))

        self.assertEqual(set(submodel_0.get_children(kind='submodel', __type=Reference)),
                         set(model.references[0:2]))
        self.assertEqual(set(submodel_1.get_children(kind='submodel', __type=Reference)),
                         set([model.references[0], model.references[2]]))

        self.assertEqual(set(submodel_0.get_children(kind='submodel', __type=Parameter)),
                         set(model.parameters[0:2]))
        self.assertEqual(set(submodel_1.get_children(kind='submodel', __type=Parameter)),
                         set([model.parameters[0], model.parameters[2]]))

        self.assertEqual(submodel_0.get_children(kind='submodel', id='rxn_01'),
                         model.reactions[1:2])
        self.assertEqual(submodel_1.get_children(kind='submodel', id='model'),
                         [model])

        self.assertEqual(submodel_0.get_children(kind='submodel', id='p_0'),
                         model.parameters[0:1])
        self.assertEqual(submodel_1.get_children(kind='submodel', id='p_2'),
                         model.parameters[2:3])


class CutTestCase(unittest.TestCase):
    def test(self):
        model, _, _ = gen_model()
        submodels = model.submodels.cut(kind='submodel')
        self.assertTrue(submodels[0].is_equal(gen_submodel_0()))
        self.assertTrue(submodels[1].is_equal(gen_submodel_1()))
