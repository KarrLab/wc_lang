""" Test clipping the flux bounds for the reactions in dFBA submodels

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-29
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from wc_lang import Model, SubmodelAlgorithm, FluxUnit
from wc_lang.transform.set_finite_dfba_flux_bounds import SetFiniteDfbaFluxBoundsTransform
import mock
import unittest


class SetFiniteDfbaFluxBoundsTransformTestCase(unittest.TestCase):
    def test(self):
        model = Model()
        submodel = model.submodels.create(algorithm=SubmodelAlgorithm.dfba)
        rxn_1 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_min=None, flux_max=None)
        rxn_2 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_min=None, flux_max=None)
        rxn_3 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_min=float('nan'), flux_max=float('nan'))
        rxn_4 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_min=float('nan'), flux_max=float('nan'))
        rxn_5 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_min=-1e3, flux_max=1e3, flux_units=FluxUnit['mol g^-1 s^-1'])
        rxn_6 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_min=-1e3, flux_max=1e3, flux_units=FluxUnit['mol g^-1 s^-1'])
        rxn_7 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_min=-1e1, flux_max=1e1, flux_units=FluxUnit['mol g^-1 s^-1'])
        rxn_8 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_min=-1e1, flux_max=1e1, flux_units=FluxUnit['mol g^-1 s^-1'])

        transform = SetFiniteDfbaFluxBoundsTransform()
        with mock.patch('wc_lang.transform.set_finite_dfba_flux_bounds.config', {
            'dfba': {
                'flux_min_bound_reversible': -2e2,
                'flux_min_bound_irreversible': 0.,
                'flux_max_bound': 1e2,
            },
        }):
            transform.run(model)

        self.assertEqual(rxn_1.flux_min, -2e2)
        self.assertEqual(rxn_1.flux_max, 1e2)
        self.assertEqual(rxn_2.flux_min, 0)
        self.assertEqual(rxn_2.flux_max, 1e2)
        self.assertEqual(rxn_3.flux_min, -2e2)
        self.assertEqual(rxn_3.flux_max, 1e2)
        self.assertEqual(rxn_4.flux_min, 0)
        self.assertEqual(rxn_4.flux_max, 1e2)
        self.assertEqual(rxn_5.flux_min, -2e2)
        self.assertEqual(rxn_5.flux_max, 1e2)
        self.assertEqual(rxn_6.flux_min, 0)
        self.assertEqual(rxn_6.flux_max, 1e2)
        self.assertEqual(rxn_7.flux_min, -1e1)
        self.assertEqual(rxn_7.flux_max, 1e1)
        self.assertEqual(rxn_8.flux_min, 0)
        self.assertEqual(rxn_8.flux_max, 1e1)
