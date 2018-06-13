'''
:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2017-10-23
:Copyright: 2016-2017, Karr Lab
:License: MIT
'''

import unittest
import os
import re
import tokenize
import token
from io import BytesIO

from wc_lang.io import Reader
from wc_lang import (RateLawEquation, RateLaw, Reaction, Submodel, Species, Function, StopCondition,
    ObjectiveFunction, Observable, Parameter, BiomassReaction, Compartment)
from wc_lang.rate_law_utils import RateLawUtils, ExpressionUtils, TokCodes, WcLangToken


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
        parameters = {}
        for specie in self.model.get_species():
            try:
                concentrations[specie.serialize()] = specie.concentration.value
            except:
                pass
        for parameter in self.model.get_parameters():
            try:
                parameters[parameter.id] = parameter.value
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
            rates = RateLawUtils.eval_reaction_rate_laws(reaction, concentrations, parameters)
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
            RateLawUtils.eval_reaction_rate_laws(reaction, {}, {})
        rate_law_equation.transcoded = 'cos(1.)'
        with self.assertRaises(NameError):
            RateLawUtils.eval_reaction_rate_laws(reaction, {}, {})
        rate_law_equation.transcoded = 'log(1.)'
        self.assertEqual(RateLawUtils.eval_reaction_rate_laws(reaction, {}, {}), [0])

        with self.assertRaisesRegexp(Exception, 'Error: unable to eval transcoded rate law'):
            RateLawUtils.eval_rate_law(RateLaw(), {'x': 1.}, {}, transcoded_equation='"x" + concentrations["x"]')


class TestExpressionUtils(unittest.TestCase):

    @staticmethod
    def get_tokens(expr):
        try:
            tokens = list(tokenize.tokenize(BytesIO(expr.encode('utf-8')).readline))
        except tokenize.TokenError:
            return []
        # strip the leading ENCODING and trailing ENDMARKER tokens
        return tokens[1:-1]

    def test_match_tokens(self):
        match_toks = ExpressionUtils.match_tokens
        single_name_pattern = (token.NAME, )
        self.assertFalse(match_toks([], []))
        self.assertFalse(match_toks(single_name_pattern, []))
        self.assertTrue(match_toks(single_name_pattern, self.get_tokens('ID2')))
        self.assertTrue(match_toks(single_name_pattern, self.get_tokens('ID3 5')))
        species_pattern = Species.Meta.token_pattern
        self.assertTrue(match_toks(species_pattern, self.get_tokens('sp1[c1]')))
        self.assertTrue(match_toks(species_pattern, self.get_tokens('sp1[c1] hi 7')))
        # note that whitespace is allowed between tokens
        self.assertTrue(match_toks(species_pattern, self.get_tokens('sp1 [ c1 ] ')))
        self.assertFalse(match_toks(species_pattern, self.get_tokens('sp1[c1')))

    def test_bad_tokens(self):
        _, errors = ExpressionUtils.deserialize(Species, 'test', '+= *= @= : {}', {})
        for bad_tok in ['+=', '*=', '@=', ':', '{', '}']:
            self.assertRegex(errors[0], ".*contains bad token\(s\):.*" + re.escape(bad_tok) + ".*")
        _, errors = ExpressionUtils.deserialize(Species, 'test', ' : \n', {})
        error = errors[0]
        self.assertIn('INDENT', error)
        self.assertIn('DEDENT', error)
        self.assertIn('\n', error)

    def test_deserialize(self):
        deserialize = ExpressionUtils.deserialize
        species = {}
        for id in ['specie_2[c]', 'specie_4[c]', 'specie_5[c]', 'specie_6[c]']:
            species[id] = Species()
        parameters = {}
        for id in ['k_cat', 'k_m', 'specie_2']:
            parameters[id] = Parameter(id=id)
        objects = {
            Species: species,
            Parameter: parameters
        }
        # deserialize(model_class, attribute, expression, objects)
        # WcLangToken(tok_code, token_string, model_type)
        expected_wc_tokens = [WcLangToken(TokCodes.other, '3'),
            WcLangToken(TokCodes.other, '*'),
            WcLangToken(TokCodes.wc_lang_obj_id, 'specie_5[c]', Species)]

        expected_modifiers = {
            Species: ['specie_5[c]'],
            Parameter: []
        }
        wc_expr_tokens, modifiers = deserialize(RateLawEquation, 'expression', '3 * specie_5[c]', objects)
        self.assertEqual(wc_expr_tokens, expected_wc_tokens)
        print('modifiers', modifiers)
        self.assertEqual(modifiers, expected_modifiers)

        wc_expr_tokens, modifiers = deserialize(RateLawEquation, 'expression', '2 * x[c]', objects)
        print(wc_expr_tokens, modifiers)
        wc_expr_tokens, modifiers = deserialize(RateLawEquation, 'expression', 'specie_2[c]', objects)
        print(wc_expr_tokens, modifiers)
        objects = {
            Species: species,
        }
        wc_expr_tokens, modifiers = deserialize(RateLawEquation, 'expression', '2 * x', objects)
        print(wc_expr_tokens, modifiers)
        '''
        expressions = ['k_cat * specie_2[c]', 'specie_4[c] / k_m', 'specie_5[c] + specie_6[c]']
        for e in expressions:
            print(deserialize(RateLawEquation, 'expression', e, objects))
        '''

    def test_eval_expr(self):
        pass
