'''
:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2017-10-23
:Copyright: 2016-2017, Karr Lab
:License: MIT
'''

import unittest
import os

from wc_lang.io import Reader
from wc_lang.core import RateLawEquation, RateLaw, Reaction, Submodel, Species
from wc_lang.rate_law_utils import RateLawUtils


class TestRateLawUtils(unittest.TestCase):

    # test_model_bad_species_names.xlsx contains the species names 'specie_1' and 'xspecie_1'.
    # the former is a prefix of the latter and would fail to be transcoded by the RE method
    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures',
                                  'test_model_bad_species_names.xlsx')

    def setUp(self):
        self.model = Reader().run(self.MODEL_FILENAME)

    def test_transcode_and_eval_rate_laws(self):

        # transcode rate laws
        RateLawUtils.transcode_rate_laws(self.model)
        concentrations = {}
        for specie in self.model.get_species():
            try:
                concentrations[specie.serialize()] = specie.concentration.value
            except:
                pass

        # evaluate the rate laws
        expected = {}
        expected['reaction_1'] = [0.0002]
        expected['reaction_2'] = [1.]
        expected['reaction_3'] = [.5, 0.003]
        expected['reaction_4'] = [0.0005]
        expected['biomass'] = []
        for reaction in self.model.get_reactions():
            rates = RateLawUtils.eval_reaction_rate_laws(reaction, concentrations)
            self.assertEqual(rates, expected[reaction.id])

    def test_eval_rate_law_exceptions(self):
        rate_law_equation = RateLawEquation(
            expression='',
            transcoded='',
        )
        rate_law = RateLaw(
            equation=rate_law_equation,
        )
        rate_law_equation.rate_law = rate_law
        reaction = Reaction(
            id='test_reaction',
            name='test_reaction',
            rate_laws=[rate_law]
        )
        rate_law_equation.transcoded = 'foo foo'
        with self.assertRaises(ValueError):
            RateLawUtils.eval_reaction_rate_laws(reaction, {})
        rate_law_equation.transcoded = 'cos(1.)'
        with self.assertRaises(NameError):
            RateLawUtils.eval_reaction_rate_laws(reaction, {})
        rate_law_equation.transcoded = 'log(1.)'
        self.assertEqual(RateLawUtils.eval_reaction_rate_laws(reaction, {}), [0])

        with self.assertRaisesRegexp(Exception, 'Error: unable to eval transcoded rate law'):
            RateLawUtils.eval_rate_law(RateLaw(), {'x': 1.}, transcoded_equation='"x" + concentrations["x"]')
