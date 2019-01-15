""" Test clipping the flux bounds for the reactions in dFBA submodels

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-29
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_lang import Model, ReactionFluxBoundUnit
from wc_lang.transform.set_finite_dfba_flux_bounds import SetFiniteDfbaFluxBoundsTransform
from wc_utils.util.ontology import wcm_ontology
import mock
import unittest


class SetFiniteDfbaFluxBoundsTransformTestCase(unittest.TestCase):
    def test(self):
        model = Model()
        submodel = model.submodels.create(framework=wcm_ontology['WCM:dynamic_flux_balance_analysis'])
        rxn_1 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_min=None, flux_max=None)
        rxn_2 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_min=None, flux_max=None)
        rxn_3 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_min=float('nan'), flux_max=float('nan'))
        rxn_4 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_min=float('nan'), flux_max=float('nan'))
        rxn_5 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_min=-1e3, flux_max=1e3, flux_bound_units=ReactionFluxBoundUnit['M s^-1'])
        rxn_6 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_min=-1e3, flux_max=1e3, flux_bound_units=ReactionFluxBoundUnit['M s^-1'])
        rxn_7 = model.reactions.create(submodel=submodel, reversible=True,
                                       flux_min=-1e1, flux_max=1e1, flux_bound_units=ReactionFluxBoundUnit['M s^-1'])
        rxn_8 = model.reactions.create(submodel=submodel, reversible=False,
                                       flux_min=-1e1, flux_max=1e1, flux_bound_units=ReactionFluxBoundUnit['M s^-1'])

        transform = SetFiniteDfbaFluxBoundsTransform()
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_lang__DOT__dfba__DOT__flux_min_bound_reversible', '-2e2')
        env.set('CONFIG__DOT__wc_lang__DOT__dfba__DOT__flux_min_bound_irreversible', '0.')
        env.set('CONFIG__DOT__wc_lang__DOT__dfba__DOT__flux_max_bound', '1e2')
        with env:
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
