""" Test creation of implicit initial distributions of zero concentrations in a model

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-29
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from wc_lang import Model, Species
from wc_lang.transform.create_implicit_dist_zero_init_concs \
    import CreateImplicitDistributionZeroInitConcentrationsTransform
import unittest


class CreateImplicitDistributionZeroInitConcentrationsTransformTestCase(unittest.TestCase):
    def test(self):
        model = Model()
        c_1 = model.compartments.create(id='c_1')
        c_2 = model.compartments.create(id='c_2')
        st_1 = model.species_types.create(id='st_1')
        st_2 = model.species_types.create(id='st_2')
        st_1_c_1 = model.species.create(species_type=st_1, compartment=c_1)
        st_1_c_2 = model.species.create(species_type=st_1, compartment=c_2)
        st_2_c_1 = model.species.create(species_type=st_2, compartment=c_1)
        st_2_c_2 = model.species.create(species_type=st_2, compartment=c_2)
        st_1_c_1.id = st_1_c_1.gen_id()
        st_1_c_2.id = st_1_c_2.gen_id()
        st_2_c_1.id = st_2_c_1.gen_id()
        st_2_c_2.id = st_2_c_2.gen_id()
        model.distribution_init_concentrations.create(species=st_1_c_1, mean=1., std=2.)
        model.distribution_init_concentrations.create(species=st_2_c_2, mean=1., std=2.)

        transform = CreateImplicitDistributionZeroInitConcentrationsTransform()
        transform.run(model)

        self.assertEqual(len(model.distribution_init_concentrations), 4)
        self.assertEqual(st_1_c_1.distribution_init_concentration.mean, 1.)
        self.assertEqual(st_1_c_2.distribution_init_concentration.mean, 0.)
        self.assertEqual(st_2_c_1.distribution_init_concentration.mean, 0.)
        self.assertEqual(st_2_c_2.distribution_init_concentration.mean, 1.)
        self.assertEqual(st_1_c_1.distribution_init_concentration.std, 2.)
        self.assertEqual(st_1_c_2.distribution_init_concentration.std, 0.)
        self.assertEqual(st_2_c_1.distribution_init_concentration.std, 0.)
        self.assertEqual(st_2_c_2.distribution_init_concentration.std, 2.)
