'''
:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2017-10-23
:Copyright: 2016-2017, Karr Lab
:License: MIT
'''

import unittest
from unittest.mock import MagicMock
import os
import re
import tokenize
import token
from io import BytesIO


from wc_lang.io import Reader
from wc_lang import (RateLawEquation, RateLaw, Reaction, Submodel, SpeciesType, Species, Function,
    StopCondition, ObjectiveFunction, Observable, Parameter, BiomassReaction, Compartment, RateLawEquation)
from wc_lang.expression_utils import (RateLawUtils, TokCodes, WcLangToken, LexMatch,
    WcLangExpression, WcLangExpressionError)


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


class TestWcLangExpression(unittest.TestCase):

    def setUp(self):
        self.objects = {
            Species: {'test_id[c]':Species(), 'x_id[c]':Species()},
            Parameter: {'test_id':Parameter(), 'param_id':Parameter()},
            Observable: {'test_id':Observable(), 'obs_id':Observable()},
            Function: {'fun_1':Function(), 'fun_2':Function()}
        }
        # also try
        self.objects_hard = {
            Species: {'test_id[c]':Species(), 'x_id[c]':Species()},
            Parameter: {'Observable':Parameter(), 'duped_id':Parameter()},
            Observable: {'test_id':Observable(), 'duped_id':Observable()},
            Function: {'Observable':Function(), 'fun_2':Function()}
        }

    @staticmethod
    def esc_re_center(re_list):
        return '.*' + '.*'.join([re.escape(an_re) for an_re in re_list]) + '.*'

    def make_wc_lang_expr(self, expr, obj_type=Species):
        return WcLangExpression(obj_type, 'expr_attr', expr, self.objects)

    def test_wc_lang_expression(self):
        expr = '3 + 5 * 6'
        wc_lang_expr = WcLangExpression(RateLawEquation, 'attr', ' ' + expr + ' ', self.objects)
        self.assertEqual(expr, wc_lang_expr.expression)
        n = 5
        wc_lang_expr = WcLangExpression(RateLawEquation, 'attr', ' + ' * n, self.objects)
        self.assertEqual([token.PLUS] * n, [tok.exact_type for tok in wc_lang_expr.tokens])
        wc_lang_expr = WcLangExpression(RateLawEquation, 'attr', '', {})
        self.assertEqual(wc_lang_expr.valid_functions, set())
        wc_lang_expr = WcLangExpression(RateLawEquation, 'attr', '', {Function:{}})
        self.assertEqual(wc_lang_expr.valid_functions, set(Function.Meta.valid_functions))
        wc_lang_expr = WcLangExpression(RateLawEquation, 'attr', '', {Function:{},ObjectiveFunction:{}})
        self.assertEqual(wc_lang_expr.valid_functions, set(Function.Meta.valid_functions))
        expr = 'id1[id2'
        with self.assertRaisesRegexp(WcLangExpressionError,
            "parsing '{}'.*creates a Python syntax error.*".format(re.escape(expr))):
            self.make_wc_lang_expr(expr)

    def test_get_wc_lang_model_type(self):
        wc_lang_expr = WcLangExpression(RateLawEquation, None, 'expr', self.objects)
        self.assertEqual(None, wc_lang_expr.get_wc_lang_model_type('NoSuchType'))
        self.assertEqual(Parameter, wc_lang_expr.get_wc_lang_model_type('Parameter'))
        self.assertEqual(Observable, wc_lang_expr.get_wc_lang_model_type('Observable'))

    def do_match_tokens_test(self, expr, pattern, expected, idx=0):
        wc_lang_expr = self.make_wc_lang_expr(expr)
        self.assertEqual(wc_lang_expr.match_tokens(pattern, idx), expected)

    def test_match_tokens(self):
        self.do_match_tokens_test('', [], False)
        single_name_pattern = (token.NAME, )
        self.do_match_tokens_test('', single_name_pattern, False)
        self.do_match_tokens_test('ID2', single_name_pattern, 'ID2')
        self.do_match_tokens_test('ID3 5', single_name_pattern, 'ID3')
        # fail to match tokens
        self.do_match_tokens_test('+ 5', single_name_pattern, False)
        # call match_tokens with 0<idx
        self.do_match_tokens_test('7 ID3', single_name_pattern, 'ID3', idx=1)
        self.do_match_tokens_test('2+ 5', single_name_pattern, False, idx=1)

        species_pattern = Species.Meta.token_pattern
        self.do_match_tokens_test('sp1[c1]+', species_pattern, 'sp1[c1]')
        self.do_match_tokens_test('sp1 +', species_pattern, False)
        # whitespace is not allowed between tokens in an ID
        self.do_match_tokens_test('sp1 [ c1 ] ', species_pattern, False)

    def do_disambiguated_id_error_test(self, expr, expected):
        wc_lang_expr = self.make_wc_lang_expr(expr)
        result = wc_lang_expr.disambiguated_id(0)
        self.assertTrue(isinstance(result, str))
        self.assertIn(expected.format(expr), result)

    def do_disambiguated_id_test(self, expr, disambig_type, id, pattern):
        wc_lang_expr = self.make_wc_lang_expr(expr)
        lex_match = wc_lang_expr.disambiguated_id(0)
        self.assertTrue(isinstance(lex_match, LexMatch))
        self.assertEqual(lex_match.num_py_tokens, len(pattern))
        self.assertEqual(len(lex_match.wc_lang_tokens), 1)
        wc_lang_token = lex_match.wc_lang_tokens[0]
        self.assertEqual(wc_lang_token,
            WcLangToken(TokCodes.wc_lang_obj_id, expr, disambig_type, id))

    def test_disambiguated_id(self):
        self.do_disambiguated_id_error_test('NotFunction.foo()',
            "contains '{}', which doesn't use 'Function' as a disambiguation model type")
        self.do_disambiguated_id_error_test('Function.foo2()',
            "contains '{}', which doesn't refer to a Function in 'objects'")

        self.do_disambiguated_id_test('Function.fun_1()', Function, 'fun_1',
            WcLangExpression.fun_type_disambig_patttern)

        self.do_disambiguated_id_error_test('NotFunction.foo()',
            "contains '{}', which doesn't use 'Function' as a disambiguation model type")
        self.do_disambiguated_id_error_test('Function.fun_1',
            "contains '{}', which uses 'Function' as a disambiguation model type but doesn't use Function syntax")
        self.do_disambiguated_id_error_test('NoSuchModel.fun_1',
            "contains '{}', but the disambiguation model type 'NoSuchModel' cannot be referenced by "
                "'Species' expressions")
        self.do_disambiguated_id_error_test('Parameter.fun_1',
            "contains '{}', but 'fun_1' is not the id of a 'Parameter'")

        self.do_disambiguated_id_test('Observable.test_id', Observable, 'test_id',
            WcLangExpression.model_type_disambig_pattern)

        # do not find a match
        wc_lang_expr = self.make_wc_lang_expr('3 * 2')
        self.assertEqual(wc_lang_expr.disambiguated_id(0), None)

    def do_related_object_id_error_test(self, expr, expected_error):
        wc_lang_expr = self.make_wc_lang_expr(expr)
        result = wc_lang_expr.related_object_id(0)
        self.assertTrue(isinstance(result, str))
        self.assertRegex(result, self.esc_re_center(expected_error))

    def test_related_object_id_errors(self):
        self.do_related_object_id_error_test('x[c]',
            ["contains the identifier(s)", "which aren't the id(s) of an object"])

    def test_related_object_id_mult_matches_error(self):
        del self.objects[Species]
        self.do_related_object_id_error_test('test_id',
            ["multiple model object id matches: 'test_id' as a Observable id, 'test_id' as a Parameter id"])

    def do_related_object_id_test(self, expr, expected_related_type, expected_id, pattern):
        wc_lang_expr = self.make_wc_lang_expr(expr)
        lex_match = wc_lang_expr.related_object_id(0)
        self.assertTrue(isinstance(lex_match, LexMatch))
        self.assertEqual(lex_match.num_py_tokens, len(pattern))
        self.assertEqual(len(lex_match.wc_lang_tokens), 1)
        wc_lang_token = lex_match.wc_lang_tokens[0]
        self.assertEqual(wc_lang_token,
            WcLangToken(TokCodes.wc_lang_obj_id, expected_id, expected_related_type, expected_id))

    def test_related_object_id_matches(self):
        self.do_related_object_id_test('test_id[c] + 3*x', Species, 'test_id[c]', Species.Meta.token_pattern)
        self.do_related_object_id_test('param_id', Parameter, 'param_id', (token.NAME, ))

        # no token matches
        wc_lang_expr = self.make_wc_lang_expr("3 * 4")
        self.assertEqual(wc_lang_expr.related_object_id(0), None)

    def do_fun_call_error_test(self, expr, expected_error, obj_type=Species):
        wc_lang_expr = self.make_wc_lang_expr(expr, obj_type=obj_type)
        result = wc_lang_expr.fun_call_id(0)
        self.assertTrue(isinstance(result, str))
        self.assertRegex(result, self.esc_re_center(expected_error))

    def test_fun_call_id_errors(self):
        self.do_fun_call_error_test('log(3)', ["contains the func name '",
            "but {}.Meta doesn't define 'valid_functions'".format(Species.__name__)],
            obj_type=Species)
        self.do_fun_call_error_test('foo(3)', ["contains the func name ",
            "but it isn't in {}.Meta.valid_functions".format(RateLawEquation.__name__)],
            obj_type=RateLawEquation)

    def test_fun_call_id(self):
        wc_lang_expr = self.make_wc_lang_expr('log(3)', obj_type=RateLawEquation)
        lex_match = wc_lang_expr.fun_call_id(0)
        self.assertTrue(isinstance(lex_match, LexMatch))
        self.assertEqual(lex_match.num_py_tokens, len(wc_lang_expr.function_pattern))
        self.assertEqual(len(lex_match.wc_lang_tokens), 2)
        self.assertEqual(lex_match.wc_lang_tokens[0], WcLangToken(TokCodes.math_fun_id, 'log'))
        self.assertEqual(lex_match.wc_lang_tokens[1], WcLangToken(TokCodes.other, '('))

        # no token match
        wc_lang_expr = self.make_wc_lang_expr('no_fun + 3')
        self.assertEqual(wc_lang_expr.fun_call_id(0), None)

    def test_bad_tokens(self):
        rv, errors = WcLangExpression(Species, 'test', '+= *= @= : {}', {}).deserialize()
        self.assertEqual(rv, None)
        for bad_tok in ['+=', '*=', '@=', ':', '{', '}']:
            self.assertRegex(errors[0], ".*contains bad token\(s\):.*" + re.escape(bad_tok) + ".*")
        # test bad tokens that don't have string values
        rv, errors = WcLangExpression(Species, 'test', """
 3
 +1""", {}).deserialize()
        self.assertEqual(rv, None)
        self.assertRegex(errors[0], re.escape("contains bad token(s)"))

    def do_deserialize_id_test(self, expr, expected_wc_tokens, expected_related_objs, model_type=RateLawEquation,
        test_objects=None):
        if test_objects is None:
            test_objects = self.objects_hard
        wc_lang_expr = WcLangExpression(model_type, 'attr', expr, test_objects)
        wc_tokens, related_objects = wc_lang_expr.deserialize()
        self.assertEqual(wc_lang_expr.errors, [])
        for obj_types in test_objects:
            if obj_types in expected_related_objs.keys():
                self.assertEqual(related_objects[obj_types], expected_related_objs[obj_types])
            else:
                self.assertEqual(related_objects[obj_types], [])

    def test_non_identifier_tokens(self):
        for expr in ['3', ' 7 * ( 5 - 3 ) / 2']:
            expected_wc_tokens = [WcLangToken(tok_code=TokCodes.other, token_string=tok) for tok in expr.strip().split()]
            self.do_deserialize_id_test(expr, expected_wc_tokens, {})

    def test_deserialize_w_ids(self):
        # test related_object_id
        expr = 'test_id'
        expected_wc_tokens = \
            [WcLangToken(tok_code=TokCodes.wc_lang_obj_id, token_string=expr, model_type=Observable, model_id=expr)]
        expected_related_objs = {Observable:[expr]}
        self.do_deserialize_id_test(expr, expected_wc_tokens, expected_related_objs)
        # test disambiguated_id
        expr = 'Parameter.duped_id + 2*Observable.duped_id'
        expected_wc_tokens = [
            WcLangToken(tok_code=TokCodes.wc_lang_obj_id, token_string='Parameter.duped_id', model_type=Parameter,
                model_id='duped_id'),
            WcLangToken(TokCodes.other, '+'),
            WcLangToken(TokCodes.other, '2'),
            WcLangToken(TokCodes.other, '*'),
            WcLangToken(tok_code=TokCodes.wc_lang_obj_id, token_string='Observable.duped_id', model_type=Observable,
                model_id='duped_id'),
        ]
        expected_related_objs = {
            Parameter:['duped_id'],
            Observable:['duped_id'],}
        self.do_deserialize_id_test(expr, expected_wc_tokens, expected_related_objs)
        # test fun_call_id
        expr = 'log(3) + fun_2() - Function.Observable()'
        expected_wc_tokens = [
            WcLangToken(tok_code=TokCodes.math_fun_id, token_string='log'),
            WcLangToken(TokCodes.other, '('),
            WcLangToken(TokCodes.other, '3'),
            WcLangToken(TokCodes.other, ')'),
            WcLangToken(TokCodes.other, '+'),
            WcLangToken(tok_code=TokCodes.wc_lang_obj_id, token_string='fun_2()', model_type=Function, model_id='fun_2'),
            WcLangToken(TokCodes.other, '-'),
            WcLangToken(tok_code=TokCodes.wc_lang_obj_id, token_string='Function.Observable()', model_type=Function,
                model_id='Observable')
        ]
        expected_related_objs = {
            Function:['fun_2', 'Observable']}
        self.do_deserialize_id_test(expr, expected_wc_tokens, expected_related_objs)

    def test_deserialize_w_multiple_ids(self):
        # at idx==0 match more than one of these related_object_id(), disambiguated_id(), fun_call_id()
        # test related_object_id and disambiguated_id'
        test_objects = {
            Parameter: {'Observable':Parameter()},
            Observable: {'test_id':Observable()}
        }
        expr = 'Observable.test_id'
        expected_wc_tokens = [
            WcLangToken(tok_code=TokCodes.wc_lang_obj_id, token_string=expr, model_type=Observable, model_id='test_id')
        ]
        expected_related_objs = {
            Observable:['test_id']}
        self.do_deserialize_id_test(expr, expected_wc_tokens, expected_related_objs, test_objects=test_objects)

        # test related_object_id and fun_call_id'
        test_objects = {
            Parameter: {'Function':Parameter()},
            Function: {'fun_2':Function()}
        }
        expr = 'Function.fun_2()'
        expected_wc_tokens = [
            WcLangToken(tok_code=TokCodes.wc_lang_obj_id, token_string=expr, model_type=Function, model_id='fun_2')
        ]
        expected_related_objs = {
            Function:['fun_2']}
        self.do_deserialize_id_test(expr, expected_wc_tokens, expected_related_objs, test_objects=test_objects)

    def do_deserialize_error_test(self, expr, expected_errors, model_type=RateLawEquation, test_objects=None):
        if test_objects is None:
            test_objects = self.objects_hard
        wc_lang_expr = WcLangExpression(model_type, 'attr', expr, test_objects)
        sb_none, errors = wc_lang_expr.deserialize()
        self.assertEqual(sb_none, None)
        # expected_errors is a list of lists of strings that should match the actual errors
        expected_errors = [self.esc_re_center(ee) for ee in expected_errors]
        self.assertEqual(len(errors), len(expected_errors),
            "Counts differ: num errors {} != Num expected errors {}".format(len(errors), len(expected_errors)))
        expected_errors_found = {}
        for expected_error in expected_errors:
            expected_errors_found[expected_error] = False
        for error in errors:
            for expected_error in expected_errors:
                if re.match(expected_error, error):
                    if expected_errors_found[expected_error]:
                        self.fail("Expected error '{}' matches again".format(expected_error))
                    expected_errors_found[expected_error] = True
        for expected_error, status in expected_errors_found.items():
            self.assertTrue(status, "Expected error '{}' not found in errors".format(expected_error))

    def test_deserialize_errors(self):
        bad_id = 'no_such_id'
        self.do_deserialize_error_test(bad_id,
            [["contains the identifier(s) '{}', which aren't the id(s) of an object".format(bad_id)]])
        bad_id = 'Observable.no_such_observable'
        self.do_deserialize_error_test(bad_id,
            [["contains multiple model object id matches: 'Observable' as a Function id, 'Observable' as a Parameter id"],
            ["contains '{}', but '{}'".format(bad_id, bad_id.split('.')[1]), "is not the id of a"]])
        bad_id = 'no_such_function'
        bad_fn_name = bad_id+'()'
        self.do_deserialize_error_test(bad_fn_name,
            [["contains the identifier(s) '{}', which aren't the id(s) of an object".format(bad_id)],
            ["contains the func name '{}', but it isn't in ".format(bad_id), "Meta.valid_functions"]])
        bad_id = 'Function'
        bad_fn_name = bad_id+'.no_such_function2()'
        self.do_deserialize_error_test(bad_fn_name,
            [["contains the identifier(s) '{}', which aren't the id(s) of an object".format(bad_id)],
            ["contains '{}', which doesn't refer to a Function".format(bad_fn_name)]])

    def test_str(self):
        expr = 'fun_1() + Parameter.param_id'
        wc_lang_expr = self.make_wc_lang_expr(expr)
        self.assertIn(expr, str(wc_lang_expr))
        self.assertIn('errors: []', str(wc_lang_expr))
        self.assertIn('wc_tokens: []', str(wc_lang_expr))
        wc_lang_expr.deserialize()
        self.assertIn(expr, str(wc_lang_expr))
        self.assertIn('errors: []', str(wc_lang_expr))
        self.assertIn('wc_tokens: [WcLangToken', str(wc_lang_expr))

    def test_model_class_lacks_meta(self):
        class Foo(object): pass
        objects = {
            Foo: {'foo_1':Foo(), 'foo_2':Foo()}
        }
        with self.assertRaisesRegexp(WcLangExpressionError, "objects key 'Foo' is not a subclass of obj_model.Model"):
            WcLangExpression(RateLawEquation, 'expr_attr', '', objects)
        with self.assertRaisesRegexp(WcLangExpressionError, "model_class 'Foo' is not a subclass of obj_model.Model"):
            WcLangExpression(Foo, 'expr_attr', '', self.objects)

    def test_eval_expr(self):
        class MockDynamicModel(object): pass
        mock_dynamic_model = MockDynamicModel()
        related_obj_val = 3
        mock_dynamic_model.eval_dynamic_obj = MagicMock(return_value=related_obj_val)

        # test combination of TokCodes
        wc_lang_expr = self.make_wc_lang_expr('4 * param_id + pow(2, obs_id)', obj_type=RateLawEquation)
        wc_lang_expr.deserialize()
        expected_val = 4 * related_obj_val + pow(2, related_obj_val)
        evaled_val = wc_lang_expr.eval_expr(None, 0, mock_dynamic_model)
        self.assertEqual(expected_val, evaled_val)

        # test types with callable ids
        c = Compartment(id='c')
        st = SpeciesType(id='x_id')
        species = Species(compartment=c, species_type=st)
        wc_lang_expr = self.make_wc_lang_expr('x_id[c] +')
        wc_lang_expr.deserialize()
        with self.assertRaisesRegexp(WcLangExpressionError, re.escape("SyntaxError: cannot eval expression '3 +' "
            "in {} with id {}".format(Species.__name__, species.id()))):
            wc_lang_expr.eval_expr(species, 0, mock_dynamic_model)

        # test different exceptions
        # syntax error
        model_type = Parameter
        wc_lang_expr = self.make_wc_lang_expr('4 *', obj_type=model_type)
        wc_lang_expr.deserialize()
        id = 'rle_1'
        with self.assertRaisesRegexp(WcLangExpressionError, "SyntaxError: cannot eval expression .* in {} with id {}".format(
            model_type.__name__, id)):
            wc_lang_expr.eval_expr(model_type(id=id), 0, mock_dynamic_model)

        # expression that could not be serialized
        expr = 'foo(6)'
        wc_lang_expr = self.make_wc_lang_expr(expr, obj_type=model_type)
        wc_lang_expr.deserialize()
        with self.assertRaisesRegexp(WcLangExpressionError, re.escape("cannot evaluate '{}', as it not been "
            "successfully deserialized".format(expr))):
            wc_lang_expr.eval_expr(model_type(id=id), 0, mock_dynamic_model)
