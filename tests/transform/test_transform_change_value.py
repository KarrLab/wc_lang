""" Tests of change value transform.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_lang import (Model, Compartment, DistributionInitConcentration, Function, FunctionExpression, Parameter,
                     Reaction, RateLawDirection, RateLaw, RateLawExpression, ReactionFluxBoundUnit, Species)
from wc_lang.transform import ChangeValueTransform
import unittest


class ChangeValueTransformTestCase(unittest.TestCase):

    def test_compartment_density(self):
        model = Model()
        c = model.compartments.create(id='c', mean_init_volume=1.)
        e = model.compartments.create(id='e', mean_init_volume=1.)
        ChangeValueTransform((('compartments', {'id': 'c'}), 'mean_init_volume'), 2.).run(model)

        self.assertEqual(c.mean_init_volume, 2.)
        self.assertEqual(e.mean_init_volume, 1.)

    def test_function_expression(self):
        model = Model()
        f = model.functions.create(id='f', expression=FunctionExpression(expression='x'))
        g = model.functions.create(id='g', expression=FunctionExpression(expression='y'))
        ChangeValueTransform((('functions', {'id': 'g'}), 'expression', 'expression'), 'z').run(model)

        self.assertEqual(f.expression.expression, 'x')
        self.assertEqual(g.expression.expression, 'z')

    def test_parameter_units(self):
        model = Model()
        p_1 = model.parameters.create(id='p_1', units='M')
        p_2 = model.parameters.create(id='p_2', units='L')
        ChangeValueTransform((('parameters', {'id': 'p_1'}), 'units'), 'm^2').run(model)

        self.assertEqual(p_1.units, 'm^2')
        self.assertEqual(p_2.units, 'L')

    def test_parameter_value(self):
        model = Model()
        p_1 = model.parameters.create(id='p_1', value=1)
        p_2 = model.parameters.create(id='p_2', value=3)
        ChangeValueTransform((('parameters', {'id': 'p_1'}), 'value'), 2).run(model)

        self.assertEqual(p_1.value, 2)
        self.assertEqual(p_2.value, 3)

    def test_reaction_reversible(self):
        model = Model()
        s_1 = model.submodels.create(id='s_1')
        s_2 = model.submodels.create(id='s_2')
        r_1_1 = model.reactions.create(submodel=s_1, id='r_1_1', reversible=False)
        r_1_2 = model.reactions.create(submodel=s_1, id='r_1_2', reversible=False)
        r_2_1 = model.reactions.create(submodel=s_2, id='r_2_1', reversible=True)
        r_2_2 = model.reactions.create(submodel=s_2, id='r_2_2', reversible=True)
        ChangeValueTransform((('reactions', {'id': 'r_2_1'}), 'reversible'), False).run(model)

        self.assertEqual(r_1_1.reversible, False)
        self.assertEqual(r_1_2.reversible, False)
        self.assertEqual(r_2_1.reversible, False)
        self.assertEqual(r_2_2.reversible, True)

    def test_reaction_flux_min(self):
        model = Model()
        s_1 = model.submodels.create(id='s_1')
        s_2 = model.submodels.create(id='s_2')
        r_1_1 = model.reactions.create(submodel=s_1, id='r_1_1', flux_min=1, flux_bound_units=ReactionFluxBoundUnit['M s^-1'])
        r_1_2 = model.reactions.create(submodel=s_1, id='r_1_2', flux_min=2, flux_bound_units=ReactionFluxBoundUnit['M s^-1'])
        r_2_1 = model.reactions.create(submodel=s_2, id='r_2_1', flux_min=3, flux_bound_units=ReactionFluxBoundUnit['M s^-1'])
        r_2_2 = model.reactions.create(submodel=s_2, id='r_2_2', flux_min=4, flux_bound_units=ReactionFluxBoundUnit['M s^-1'])
        ChangeValueTransform((('reactions', {'id': 'r_1_2'}), 'flux_min'), 0).run(model)

        self.assertEqual(r_1_1.flux_min, 1)
        self.assertEqual(r_1_2.flux_min, 0)
        self.assertEqual(r_2_1.flux_min, 3)
        self.assertEqual(r_2_2.flux_min, 4)

    def test_reaction_flux_max(self):
        model = Model()
        s_1 = model.submodels.create(id='s_1')
        s_2 = model.submodels.create(id='s_2')
        r_1_1 = model.reactions.create(submodel=s_1, id='r_1_1', flux_max=1)
        r_1_2 = model.reactions.create(submodel=s_1, id='r_1_2', flux_max=2)
        r_2_1 = model.reactions.create(submodel=s_2, id='r_2_1', flux_max=3)
        r_2_2 = model.reactions.create(submodel=s_2, id='r_2_2', flux_max=4)
        ChangeValueTransform((('reactions', {'id': 'r_2_2'}), 'flux_max'), 0).run(model)

        self.assertEqual(r_1_1.flux_max, 1)
        self.assertEqual(r_1_2.flux_max, 2)
        self.assertEqual(r_2_1.flux_max, 3)
        self.assertEqual(r_2_2.flux_max, 0)

    def test_reaction_expression(self):
        model = Model()
        s_1 = model.submodels.create(id='s_1')
        r_1_1 = model.reactions.create(submodel=s_1, id='r_1_1')
        rl_f = model.rate_laws.create(reaction=r_1_1, direction=RateLawDirection.forward, expression=RateLawExpression(expression='x'))
        rl_b = model.rate_laws.create(reaction=r_1_1, direction=RateLawDirection.backward, expression=RateLawExpression(expression='y'))
        ChangeValueTransform((('reactions', {'id': 'r_1_1'}),
                              ('rate_laws', {'direction': RateLawDirection.forward}),
                              'expression',
                              'expression'), 'z').run(model)

        self.assertEqual(rl_f.expression.expression, 'z')
        self.assertEqual(rl_b.expression.expression, 'y')

    def test_reaction_k_cat(self):
        model = Model()
        s_1 = model.submodels.create(id='s_1')
        r_1_1 = model.reactions.create(submodel=s_1, id='r_1_1')
        k_cat_f = model.parameters.create(id='k_cat_f', value=1)
        k_cat_b = model.parameters.create(id='k_cat_b', value=2)
        rl_f = model.rate_laws.create(reaction=r_1_1, direction=RateLawDirection.forward,
                                      expression=RateLawExpression(
                                          parameters=[k_cat_f]))
        rl_b = model.rate_laws.create(reaction=r_1_1, direction=RateLawDirection.backward,
                                      expression=RateLawExpression(
                                          parameters=[k_cat_b]))
        rl_f.id = rl_f.gen_id()
        rl_b.id = rl_b.gen_id()
        ChangeValueTransform((('parameters', {'id': 'k_cat_b'}), 'value'), 0).run(model)

        self.assertEqual(k_cat_f.value, 1)
        self.assertEqual(k_cat_b.value, 0)

    def test_species_concentration_units(self):
        model = Model()
        st_1 = model.species_types.create(id='st_1')
        st_2 = model.species_types.create(id='st_2')
        c_1 = model.compartments.create(id='c_1')
        c_2 = model.compartments.create(id='c_2')
        st_1_c_1 = model.species.create(species_type=st_1, compartment=c_1)
        st_1_c_2 = model.species.create(species_type=st_1, compartment=c_2)
        st_2_c_1 = model.species.create(species_type=st_2, compartment=c_1)
        st_2_c_2 = model.species.create(species_type=st_2, compartment=c_2)
        st_1_c_1.id = st_1_c_1.gen_id()
        st_1_c_2.id = st_1_c_2.gen_id()
        st_2_c_1.id = st_2_c_1.gen_id()
        st_2_c_2.id = st_2_c_2.gen_id()

        st_1_c_1.distribution_init_concentration = DistributionInitConcentration(units='u')
        st_1_c_2.distribution_init_concentration = DistributionInitConcentration(units='v')
        st_2_c_1.distribution_init_concentration = DistributionInitConcentration(units='w')
        st_2_c_2.distribution_init_concentration = DistributionInitConcentration(units='x')
        st_1_c_1.distribution_init_concentration.id = st_1_c_1.distribution_init_concentration.gen_id()
        st_1_c_2.distribution_init_concentration.id = st_1_c_2.distribution_init_concentration.gen_id()
        st_2_c_1.distribution_init_concentration.id = st_2_c_1.distribution_init_concentration.gen_id()
        st_2_c_2.distribution_init_concentration.id = st_2_c_2.distribution_init_concentration.gen_id()

        ChangeValueTransform((('species', {'id': 'st_1[c_1]'}),
                              'distribution_init_concentration',
                              'units'), 'a').run(model)
        self.assertEqual(st_1_c_1.distribution_init_concentration.units, 'a')
        self.assertEqual(st_1_c_2.distribution_init_concentration.units, 'v')
        self.assertEqual(st_2_c_1.distribution_init_concentration.units, 'w')
        self.assertEqual(st_2_c_2.distribution_init_concentration.units, 'x')

    def test_species_concentration_value(self):
        model = Model()
        st_1 = model.species_types.create(id='st_1')
        st_2 = model.species_types.create(id='st_2')
        c_1 = model.compartments.create(id='c_1')
        c_2 = model.compartments.create(id='c_2')
        st_1_c_1 = model.species.create(species_type=st_1, compartment=c_1)
        st_1_c_2 = model.species.create(species_type=st_1, compartment=c_2)
        st_2_c_1 = model.species.create(species_type=st_2, compartment=c_1)
        st_2_c_2 = model.species.create(species_type=st_2, compartment=c_2)
        st_1_c_1.id = st_1_c_1.gen_id()
        st_1_c_2.id = st_1_c_2.gen_id()
        st_2_c_1.id = st_2_c_1.gen_id()
        st_2_c_2.id = st_2_c_2.gen_id()

        st_1_c_1.distribution_init_concentration = DistributionInitConcentration(mean=1, std=2)
        st_1_c_2.distribution_init_concentration = DistributionInitConcentration(mean=2, std=2)
        st_2_c_1.distribution_init_concentration = DistributionInitConcentration(mean=3, std=2)
        st_2_c_2.distribution_init_concentration = DistributionInitConcentration(mean=4, std=2)
        st_1_c_1.distribution_init_concentration.id = st_1_c_1.distribution_init_concentration.gen_id()
        st_1_c_2.distribution_init_concentration.id = st_1_c_2.distribution_init_concentration.gen_id()
        st_2_c_1.distribution_init_concentration.id = st_2_c_1.distribution_init_concentration.gen_id()
        st_2_c_2.distribution_init_concentration.id = st_2_c_2.distribution_init_concentration.gen_id()

        ChangeValueTransform((('species', {'id': st_2_c_1.id}),
                              'distribution_init_concentration',
                              'mean'), 0).run(model)
        self.assertEqual(st_1_c_1.distribution_init_concentration.mean, 1)
        self.assertEqual(st_1_c_2.distribution_init_concentration.mean, 2)
        self.assertEqual(st_2_c_1.distribution_init_concentration.mean, 0)
        self.assertEqual(st_2_c_2.distribution_init_concentration.mean, 4)

        ChangeValueTransform((('species', {'id': st_2_c_1.id}),
                              'distribution_init_concentration',
                              'std'), 0).run(model)
        self.assertEqual(st_1_c_1.distribution_init_concentration.std, 2)
        self.assertEqual(st_1_c_2.distribution_init_concentration.std, 2)
        self.assertEqual(st_2_c_1.distribution_init_concentration.std, 0)
        self.assertEqual(st_2_c_2.distribution_init_concentration.std, 2)

    def test_attr_path_to_from_str(self):
        t1 = ChangeValueTransform('species', 0)
        t2 = ChangeValueTransform(None, 0)
        t2.attr_path = ChangeValueTransform.attr_path_from_str(t1.attr_path_to_str())
        self.assertEqual(t1.attr_path, t2.attr_path)

        t1 = ChangeValueTransform(['species'], 0)
        t2 = ChangeValueTransform(None, 0)
        t2.attr_path = ChangeValueTransform.attr_path_from_str(t1.attr_path_to_str())
        self.assertEqual(t1.attr_path, t2.attr_path)

        t1 = ChangeValueTransform(['species', 'reactions'], 0)
        t2 = ChangeValueTransform(None, 0)
        t2.attr_path = ChangeValueTransform.attr_path_from_str(t1.attr_path_to_str())
        self.assertEqual(t1.attr_path, t2.attr_path)

        t1 = ChangeValueTransform([['species', {'id': 'st_1[c_1]'}],
                                   'distribution_init_concentration',
                                   'std'], 0)
        t2 = ChangeValueTransform(None, 0)
        t2.attr_path = ChangeValueTransform.attr_path_from_str(t1.attr_path_to_str())
        self.assertEqual(t1.attr_path, t2.attr_path)

        t1 = ChangeValueTransform(['rate_laws', {'direction': RateLawDirection.forward}], 0)
        t2 = ChangeValueTransform(None, 0)
        t2.attr_path = ChangeValueTransform.attr_path_from_str(t1.attr_path_to_str())
        self.assertEqual(t1.attr_path, t2.attr_path)

    def test___eq__(self):
        t1 = ChangeValueTransform((('species', {'id': 'st_1[c_1]'}),
                                   'distribution_init_concentration',
                                   'std'), 0)
        t2 = ChangeValueTransform((('species', {'id': 'st_1[c_1]'}),
                                   'distribution_init_concentration',
                                   'std'), 0)
        self.assertEqual(t1, t2)

        t2.value = 1.
        self.assertNotEqual(t1, t2)

        self.assertNotEqual(t1, 1.)
        self.assertNotEqual(1., t1)
