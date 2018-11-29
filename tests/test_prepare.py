""" Test WC model preparation

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-10-22
:Copyright: 2017, Karr Lab
:License: MIT
"""
from wc_lang import (Model, Submodel, Reaction, SpeciesType, Species,
                     Compartment, SubmodelAlgorithm)
from wc_lang import prepare
import mock
import os
import unittest
import wc_lang.config.core
import wc_lang.io


class PrepareModelTransformTestCase(unittest.TestCase):
    def test_run(self):
        filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'example-model.xlsx')
        model = wc_lang.io.Reader().run(filename)
        transform = prepare.PrepareModelTransform()
        transform.run(model)


class CreateImplicitZeroConcentrationsTransformTestCase(unittest.TestCase):
    def test(self):
        model = Model()
        c_1 = model.compartments.create(id='c_1')
        c_2 = model.compartments.create(id='c_2')
        st_1 = model.species_types.create(id='st_1')
        st_2 = model.species_types.create(id='st_2')
        st_1_c_1 = model.species.create(id=Species.gen_id('st_1', 'c_1'), species_type=st_1, compartment=c_1)
        st_1_c_2 = model.species.create(id=Species.gen_id('st_1', 'c_2'), species_type=st_1, compartment=c_2)
        st_2_c_1 = model.species.create(id=Species.gen_id('st_2', 'c_1'), species_type=st_2, compartment=c_1)
        st_2_c_2 = model.species.create(id=Species.gen_id('st_2', 'c_2'), species_type=st_2, compartment=c_2)
        model.concentrations.create(species=st_1_c_1, mean=1., std=2.)
        model.concentrations.create(species=st_2_c_2, mean=1., std=2.)

        transform = prepare.CreateImplicitZeroConcentrationsTransform()
        transform.run(model)

        self.assertEqual(len(model.concentrations), 4)
        self.assertEqual(st_1_c_1.concentration.mean, 1.)
        self.assertEqual(st_1_c_2.concentration.mean, 0.)
        self.assertEqual(st_2_c_1.concentration.mean, 0.)
        self.assertEqual(st_2_c_2.concentration.mean, 1.)
        self.assertEqual(st_1_c_1.concentration.std, 2.)
        self.assertEqual(st_1_c_2.concentration.std, 0.)
        self.assertEqual(st_2_c_1.concentration.std, 0.)
        self.assertEqual(st_2_c_2.concentration.std, 2.)


class CreateImplicitDfbaExchangeReactionsTransformTestCase(unittest.TestCase):
    def test(self):
        model = Model()
        submodel = model.submodels.create(id='submdl', name='submodel', algorithm=SubmodelAlgorithm.dfba)

        comps = [
            model.compartments.create(id='c', name='cytosol'),
            model.compartments.create(id='d', name='dna'),
            model.compartments.create(id='e', name='extracellular space'),
        ]
        sts = [
            model.species_types.create(id='st_1', name='species type 1'),
            model.species_types.create(id='st_2', name='species type 2'),
            model.species_types.create(id='st_3', name='species type 3'),
        ]
        specs = []
        for st in sts:
            spec_comps = []
            specs.append(spec_comps)
            for comp in comps:
                spec_comps.append(model.species.create(species_type=st, compartment=comp))

        rxn = model.reactions.create(submodel=submodel)
        rxn.participants.create(species=specs[0][0])
        rxn.participants.create(species=specs[0][1])
        rxn.participants.create(species=specs[0][2])

        rxn = model.reactions.create(submodel=submodel)
        rxn.participants.create(species=specs[1][0])
        rxn.participants.create(species=specs[1][1])
        rxn.participants.create(species=specs[1][2])

        rxn = model.reactions.create(submodel=submodel)
        rxn.participants.create(species=specs[2][0])
        rxn.participants.create(species=specs[2][1])

        transform = prepare.CreateImplicitDfbaExchangeReactionsTransform()
        transform.run(model)

        self.assertEqual(len(model.reactions), 5)
        self.assertEqual(len(submodel.reactions), 5)
        self.assertNotEqual(model.reactions.get_one(id='__dfba_ex_submdl_st_1_e'), None)
        self.assertNotEqual(model.reactions.get_one(id='__dfba_ex_submdl_st_2_e'), None)
        self.assertEqual(model.reactions.get_one(id='__dfba_ex_submdl_st_3_e'), None)

        rxn = model.reactions.get_one(id='__dfba_ex_submdl_st_1_e')
        self.assertEqual(rxn.name, 'dFBA exchange (submodel, species type 1, extracellular space)')
        self.assertEqual(len(rxn.participants), 1)
        self.assertEqual(rxn.participants[0].species, specs[0][2])
        self.assertEqual(rxn.participants[0].coefficient, 1.)
        self.assertEqual(rxn.reversible, True)


class SetFiniteDfbaFluxBoundsTransformTestCase(unittest.TestCase):
    def test(self):
        model = Model()
        submodel = model.submodels.create(algorithm=SubmodelAlgorithm.dfba)
        rxn_1 = model.reactions.create(submodel=submodel, reversible=True, min_flux=None, max_flux=None)
        rxn_2 = model.reactions.create(submodel=submodel, reversible=False, min_flux=None, max_flux=None)
        rxn_3 = model.reactions.create(submodel=submodel, reversible=True, min_flux=float('nan'), max_flux=float('nan'))
        rxn_4 = model.reactions.create(submodel=submodel, reversible=False, min_flux=float('nan'), max_flux=float('nan'))
        rxn_5 = model.reactions.create(submodel=submodel, reversible=True, min_flux=-1e3, max_flux=1e3)
        rxn_6 = model.reactions.create(submodel=submodel, reversible=False, min_flux=-1e3, max_flux=1e3)
        rxn_7 = model.reactions.create(submodel=submodel, reversible=True, min_flux=-1e1, max_flux=1e1)
        rxn_8 = model.reactions.create(submodel=submodel, reversible=False, min_flux=-1e1, max_flux=1e1)

        transform = prepare.SetFiniteDfbaFluxBoundsTransform()
        with mock.patch('wc_lang.prepare.config', {
            'dfba': {
                'min_reversible_flux_bound': -2e2,
                'min_irreversible_flux_bound': 0.,
                'max_flux_bound': 1e2,
            },
        }):
            transform.run(model)

        self.assertEqual(rxn_1.min_flux, -2e2)
        self.assertEqual(rxn_1.max_flux, 1e2)
        self.assertEqual(rxn_2.min_flux, 0)
        self.assertEqual(rxn_2.max_flux, 1e2)
        self.assertEqual(rxn_3.min_flux, -2e2)
        self.assertEqual(rxn_3.max_flux, 1e2)
        self.assertEqual(rxn_4.min_flux, 0)
        self.assertEqual(rxn_4.max_flux, 1e2)
        self.assertEqual(rxn_5.min_flux, -2e2)
        self.assertEqual(rxn_5.max_flux, 1e2)
        self.assertEqual(rxn_6.min_flux, 0)
        self.assertEqual(rxn_6.max_flux, 1e2)
        self.assertEqual(rxn_7.min_flux, -1e1)
        self.assertEqual(rxn_7.max_flux, 1e1)
        self.assertEqual(rxn_8.min_flux, 0)
        self.assertEqual(rxn_8.max_flux, 1e1)
