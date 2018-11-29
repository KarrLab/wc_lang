""" Test clipping the flux bounds for the reactions in dFBA submodels

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-29
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from wc_lang import Model, SubmodelAlgorithm
from wc_lang.transform.set_finite_dfba_flux_bounds import SetFiniteDfbaFluxBoundsTransform
import mock
import unittest


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

        transform = SetFiniteDfbaFluxBoundsTransform()
        with mock.patch('wc_lang.transform.set_finite_dfba_flux_bounds.config', {
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
