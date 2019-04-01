""" Test clipping the flux bounds for the reactions in dFBA submodels

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-29
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_lang.core import Model, FluxBounds
from wc_lang.transform.set_finite_dfba_flux_bounds import SetFiniteDfbaFluxBoundsTransform
from wc_onto import onto
from wc_utils.util.units import unit_registry
import mock
import unittest


class SetFiniteDfbaFluxBoundsTransformTestCase(unittest.TestCase):
    def test(self):
        model = Model()
        submodel = model.submodels.create(framework=onto['WC:dynamic_flux_balance_analysis'])
        rxn_1 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_bounds=None)
        rxn_2 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_bounds=FluxBounds(min=None, max=None))
        rxn_3 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_bounds=FluxBounds(min=float('nan'), max=float('nan')))
        rxn_4 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_bounds=FluxBounds(min=float('nan'), max=float('nan')))
        rxn_5 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_bounds=FluxBounds(min=-1e3, max=1e3, units=unit_registry.parse_units('M s^-1')))
        rxn_6 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_bounds=FluxBounds(min=-1e3, max=1e3, units=unit_registry.parse_units('M s^-1')))
        rxn_7 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_bounds=FluxBounds(min=-1e1, max=1e1, units=unit_registry.parse_units('M s^-1')))
        rxn_8 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_bounds=FluxBounds(min=-1e1, max=1e1, units=unit_registry.parse_units('M s^-1')))

        transform = SetFiniteDfbaFluxBoundsTransform()
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_lang__DOT__dfba__DOT__flux_bounds__DOT__min_reversible', '-2e2')
        env.set('CONFIG__DOT__wc_lang__DOT__dfba__DOT__flux_bounds__DOT__min_irreversible', '0.')
        env.set('CONFIG__DOT__wc_lang__DOT__dfba__DOT__flux_bounds__DOT__max', '1e2')
        with env:
            transform.run(model)

        self.assertEqual(rxn_1.flux_bounds.min, -2e2)
        self.assertEqual(rxn_1.flux_bounds.max, 1e2)
        self.assertEqual(rxn_2.flux_bounds.min, 0)
        self.assertEqual(rxn_2.flux_bounds.max, 1e2)
        self.assertEqual(rxn_3.flux_bounds.min, -2e2)
        self.assertEqual(rxn_3.flux_bounds.max, 1e2)
        self.assertEqual(rxn_4.flux_bounds.min, 0)
        self.assertEqual(rxn_4.flux_bounds.max, 1e2)
        self.assertEqual(rxn_5.flux_bounds.min, -2e2)
        self.assertEqual(rxn_5.flux_bounds.max, 1e2)
        self.assertEqual(rxn_6.flux_bounds.min, 0)
        self.assertEqual(rxn_6.flux_bounds.max, 1e2)
        self.assertEqual(rxn_7.flux_bounds.min, -1e1)
        self.assertEqual(rxn_7.flux_bounds.max, 1e1)
        self.assertEqual(rxn_8.flux_bounds.min, 0)
        self.assertEqual(rxn_8.flux_bounds.max, 1e1)
