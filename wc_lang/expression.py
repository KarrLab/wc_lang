""" Utilities for processing mathematical expressions used by wc_lang models

:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2017-10-23
:Copyright: 2016-2018, Karr Lab
:License: MIT
"""
import token
import tokenize
import obj_model
import wc_lang.core
from collections import namedtuple
from io import BytesIO
from wc_utils.util.enumerate import CaseInsensitiveEnum
from wc_utils.util.misc import DFSMAcceptor

# TODOS
'''
replace MakeModels
'''
'''
build
wc_lang:
    right template for wc_sim_test
    use ParsedExpression to deserialize and validate all wc_lang expressions:
        last, convert DfbaObjective, and RateLawExpression (what about k_cat and k_m?)
    fix test_io_roundtrip.py:test_create()
    in DynamicReaction() init submodels must have compartments: enforce this
    more rigerously create rate laws, whose kinetics must be MA, MM, or custom
wc_sim:
    put dict of submodels in dynamic_model
    give abstract compartments control over their volumes
    can a reaction's reactants and modifiers be contained in multiple compartements
    ensure that all related objects made before their using objects in DynamicModel.__init__()
    handle any multiplicity of StopConditions; allow execution to choose among multiple StopConditions,
        or use them all
extra:
    Parameter: make Parameter.id unique; ensure that Parameters are constants
    make test for "AssertionError: No matches or errors found in 'obs_1'" by providing empty objects
    test multiple instances of the same used model in an expression
    incorporate id into error_suffix in test_eval
    test with real WC models
    perhaps make model_class a model, so that test_eval can provide ids
    rename Function to Macro; need to change model files too -- argh
    expand Jupyter example
    support logical operators (or, and, not) in expressions, esp. StopConditions
cleanup:
    robust tests of run_results that examine the data
    more docstrings and better naming for ExpressionVerifier
    use ExpressionMethods.make_expression_obj() wherever possible
    better error than "Not a linear expression of species and observables"
    do a better job of getting the plural in "if not used_model_type_attr.endswith('s'):"
    have valid_functions defined as sets, not tuples
    replace all validation and deserialization code
    fix test_find_shared_species
'''


class TokCodes(int, CaseInsensitiveEnum):
    """ Token codes used in parsed expressions """
    wc_lang_obj_id = 1
    math_fun_id = 2
    number = 3
    op = 4
    other = 5


# a matched token pattern used by tokenize
IdMatch = namedtuple('IdMatch', 'model_type, token_pattern, match_string')
IdMatch.__doc__ += ': Matched token pattern used by tokenize'
IdMatch.model_type.__doc__ = 'The type of Model matched'
IdMatch.token_pattern.__doc__ = 'The token pattern used by the match'
IdMatch.match_string.__doc__ = 'The matched string'


# a token in a parsed wc_lang expression, returned in a list by tokenize
WcLangToken = namedtuple('WcLangToken', 'tok_code, token_string, model_type, model_id, model')
# make model_type, model_id, and model optional: see https://stackoverflow.com/a/18348004
WcLangToken.__new__.__defaults__ = (None, None, None)
WcLangToken.__doc__ += ': Token in a parsed wc_lang expression'
WcLangToken.tok_code.__doc__ = 'TokCodes encoding'
WcLangToken.token_string.__doc__ = "The token's string"
WcLangToken.model_type.__doc__ = "When tok_code is wc_lang_obj_id, the wc_lang obj's type"
WcLangToken.model_id.__doc__ = "When tok_code is wc_lang_obj_id, the wc_lang obj's id"
WcLangToken.model.__doc__ = "When tok_code is wc_lang_obj_id, the wc_lang obj"


# result returned by a tokens lexer, like disambiguated_id()
LexMatch = namedtuple('LexMatch', 'wc_lang_tokens, num_py_tokens')
LexMatch.__doc__ += ': result returned by a lexer method that matches a wc_lang expression element'
LexMatch.wc_lang_tokens.__doc__ = 'List of WcLangTokens created'
LexMatch.num_py_tokens.__doc__ = 'Number of Python tokens consumed'


class ParsedExpressionError(Exception):
    """ Exception raised for errors in `ParsedExpression`

    Attributes:
        message (:obj:`str`): the exception's message
    """

    def __init__(self, message=None):
        super().__init__(message)


class ParsedExpression(object):
    """ An expression in a wc_lang Model

    Expressions are currently (July, 2018) used in five `wc_lang` `Model`s: `RateLawExpression`, `Function`,
    `StopCondition` (which is just a special case of `Function` that returns a boolean), `DfbaObjective`,
    and `Observable`. These expressions are limited Python expressions with specific semantics:

    * They must be syntactically correct Python.
    * No Python keywords, strings, or tokens that do not belong in expressions are allowed.
    * All Python identifiers must be the ID of an object in a whole-cell model, or components of
        an ID of an object in the model, or the name of a function in the `math` package. Objects in the model
        are provided in `_objects`, and the allowed subset of functions in `math` must be provided in an
        iterator in the `valid_functions` attribute of the `Meta` class of a model whose whose expression
        is being processed.
    * Currently (July, 2018), identifiers may refer to `Species`s, `Parameter`s, `Observable`s, `Reaction`s,
        `Observable`'s and `BiomassReaction`s.
    * Cycles of references are illegal.
    * An identifier must unambiguously refer to exactly one related `Model` in a model.
    * Each `wc_lang` `Model` that can be used in an expression must have an ID that is a simple Python identifier,
        or define `token_pattern` as an attribute that describes the `Model`'s syntactic Python structure. See
        `Species` for an example.
    * Every expression must be computable at any time during a simulation. The evaluation of an expression
        always occurs at a precise simulation time, which is implied by the expression but not explicitly
        represented. E.g., a reference to a `Species` means its concentration at the time the expression is
        `eval`ed. These are the meanings of references:
        * `Species`: its current concentration
        * `Parameter`: its value, which is static
        * `Observable`: its current value, whose units depend on its definition
        * `Reaction`: its current flux
        * `BiomassReaction`: its current flux
    * The modeller is responsible for ensuring that units in expressions are internally consistent and appropriate
        for the expression's use

    Attributes:
        model_class (:obj:`obj_model.Model`): the `wc_lang` `Model` which has an expression
        attribute (:obj:`str`): the attribute name of the expression in `model_class`
        expression (:obj:`str`): the expression defined in the wc_lang Model
        tokens (:obj:`list` of :obj:`namedtuple`): a list of Python tokens generated by `tokenize.tokenize()`
        _objects (:obj:`dict`): dict of wc_lang Models that might be referenced in expression; maps
            model type to a dict mapping ids to Model instances
        valid_models (:obj:`set`): wc_lang Models that `model_class` objects are allowed to use,
            as indicated in `model_class.Meta.valid_models`, intersected with `_objects.keys()`
            might be referenced in expression; maps
        valid_functions (:obj:`set`): the union of all `valid_functions` attributes for `_objects`
        related_objects (:obj:`dict`): models that are referenced in `expression`; maps model type to
            dict that maps model id to model instance
        lin_coeffs (:obj:`dict`): linear coefficients of models that are referenced in `expression`; 
            maps model type to dict that maps models to coefficients
        errors (:obj:`list` of :obj:`str`): errors found when parsing an `expression` fails
        wc_tokens (:obj:`list` of :obj:`WcLangToken`): tokens obtained when an `expression` is successfully
            `tokenize`d; if empty, then this `ParsedExpression` cannot use `eval()`
    """

    # Function.identifier()
    fun_type_disambig_patttern = (token.NAME, token.DOT, token.NAME, token.LPAR, token.RPAR)
    # ModelType.model_id
    model_type_disambig_pattern = (token.NAME, token.DOT, token.NAME)
    function_pattern = (token.NAME, token.LPAR)

    # enumerate and detect Python tokens that are illegal in wc_lang expressions
    # TODO: consider the handful of other tokens that may also be illegal: COMMA, SEMI, TILDE, CIRCUMFLEX, and AT
    illegal_token_names = ['ENDMARKER', 'NEWLINE', 'INDENT', 'DEDENT', 'COLON', 'LBRACE', 'RBRACE',
                           'PLUSEQUAL', 'MINEQUAL', 'STAREQUAL', 'SLASHEQUAL', 'PERCENTEQUAL', 'AMPEREQUAL', 'VBAREQUAL',
                           'CIRCUMFLEXEQUAL', 'LEFTSHIFTEQUAL', 'RIGHTSHIFTEQUAL', 'DOUBLESTAREQUAL', 'DOUBLESLASHEQUAL',
                           'ATEQUAL', 'RARROW', 'ELLIPSIS', 'AWAIT', 'ASYNC', 'ERRORTOKEN', 'N_TOKENS', 'NT_OFFSET']
    illegal_tokens = set()
    for illegal_name in illegal_token_names:
        illegal_tokens.add(getattr(token, illegal_name))

    def __init__(self, model_class, attribute, expression, objects):
        """ Create an instance of ParsedExpression

        Raises:
            (:obj:`ParsedExpressionError`): if `model_class` is not a subclass of `obj_model.Model`,
                or lexical analysis of `expression` raises an exception,
                or `objects` includes model types that `model_class` should not reference
        """
        if not issubclass(model_class, obj_model.Model):
            raise ParsedExpressionError("model_class '{}' is not a subclass of obj_model.Model".format(
                model_class.__name__))
        if not hasattr(model_class.Meta, 'valid_models'):
            raise ParsedExpressionError("model_class '{}' doesn't have a 'Meta.valid_models' attribute".format(
                model_class.__name__))
        self.valid_models = set()
        for valid_model_type_name in model_class.Meta.valid_models:
            valid_model_type = getattr(wc_lang.core, valid_model_type_name)
            if valid_model_type in objects:
                self.valid_models.add(valid_model_type)
        for obj_type in self.valid_models:
            if not issubclass(obj_type, obj_model.Model):   # pragma    no cover
                raise ParsedExpressionError("objects entry '{}' is not a subclass of obj_model.Model".format(
                    obj_type.__name__))
        self.valid_functions = set()
        if hasattr(model_class.Meta, 'valid_functions'):
            self.valid_functions.update(model_class.Meta.valid_functions)

        self._objects = objects
        self.model_class = model_class
        self.attribute = attribute
        # strip leading and trailing whitespace from expression, which would create a bad token error
        self.expression = expression.strip()

        try:
            g = tokenize.tokenize(BytesIO(self.expression.encode('utf-8')).readline)
            # strip the leading ENCODING and trailing ENDMARKER tokens
            self.tokens = list(g)[1:-1]
        except tokenize.TokenError as e:
            raise ParsedExpressionError("parsing '{}', a {}.{}, creates a Python syntax error: '{}'".format(
                self.expression, self.model_class.__name__, self.attribute, str(e)))

        self.__reset_tokenization()

    def __reset_tokenization(self):
        """ Reset tokenization
        """
        self.related_objects = {}
        self.lin_coeffs = {}
        for model_type in self.valid_models:
            self.related_objects[model_type] = {}
            self.lin_coeffs[model_type] = {}

        self.errors = []
        self.wc_tokens = []

    def get_wc_lang_model_type(self, model_type_name):
        """ Find the `wc_lang` model type corresponding to `model_type_name`

        Args:
            model_type_name (:obj:`str`): the name of a purported `wc_lang` model type in an expression

        Returns:
            :obj:`object`: `None` if no model named `model_type_name` exists in `self.valid_models`,
                else the type of the model with that name
        """
        for model_type in self.valid_models:
            if model_type_name == model_type.__name__:
                return model_type
        return None

    def match_tokens(self, token_pattern, idx):
        """ Indicate whether `tokens` begins with a pattern of tokens that match `token_pattern`

        Args:
            token_pattern (:obj:`tuple` of :obj:`int`): a tuple of Python token numbers, taken from the
            `token` module
            idx (:obj:`int`): current index into `tokens`

            Returns:
                :obj:`object`: :obj:`bool`, False if the initial elements of `tokens` do not match the
                syntax in `token_pattern`, or :obj:`str`, the matching string
        """
        if not token_pattern:
            return False
        if len(self.tokens)-idx < len(token_pattern):
            return False
        for tok_idx, token_pat_num in enumerate(token_pattern):
            if self.tokens[idx+tok_idx].exact_type != token_pat_num:
                return False
            # because a wc_lang ID shouldn't contain white space, do not allow it between the self.tokens
            # that match token_pattern
            if 0 < tok_idx and self.tokens[idx+tok_idx-1].end != self.tokens[idx+tok_idx].start:
                return False
        match_val = ''.join([self.tokens[idx+i].string for i in range(len(token_pattern))])
        return match_val

    def disambiguated_id(self, idx, case_fold_match=False):
        """ Try to parse a disambuated `wc_lang` id from `self.tokens` at `idx`

        Look for a disambugated id (either a Function written as `Function.identifier()`, or a
        Model written as `ModelType.model_id`). If tokens do not match, return `None`. If tokens match,
        but their values are wrong, return an error `str`.
        If a disambugated id is found, return a `LexMatch` describing it.

        Args:
            idx (:obj:`int`): current index into `tokens`
            case_fold_match (:obj:`bool`, optional): if set, `casefold()` identifiers before matching;
                in a `WcLangToken`, `token_string` retains the original expression text, while `model_id`
                contains the casefold'ed value; identifier keys in `self._objects` must already be casefold'ed;
                default=False

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a disambugated id is found, return a `LexMatch` describing it.
        """
        fun_match = self.match_tokens(self.fun_type_disambig_patttern, idx)
        if fun_match:
            possible_macro_id = self.tokens[idx+2].string
            if case_fold_match:
                possible_macro_id = possible_macro_id.casefold()
            # the disambiguation model type must be Function
            if self.tokens[idx].string != wc_lang.core.Function.__name__:
                return ("'{}', a {}.{}, contains '{}', which doesn't use 'Function' as a disambiguation "
                        "model type".format(self.expression, self.model_class.__name__, self.attribute, fun_match))
            # the identifier must be in the Function objects
            if wc_lang.core.Function not in self.valid_models or possible_macro_id not in self._objects[wc_lang.core.Function]:
                return "'{}', a {}.{}, contains '{}', which doesn't refer to a Function".format(
                    self.expression, self.model_class.__name__, self.attribute, fun_match)
            return LexMatch([WcLangToken(TokCodes.wc_lang_obj_id, fun_match, wc_lang.core.Function,
                                         possible_macro_id, self._objects[wc_lang.core.Function][possible_macro_id])],
                            len(self.fun_type_disambig_patttern))

        disambig_model_match = self.match_tokens(self.model_type_disambig_pattern, idx)
        if disambig_model_match:
            disambig_model_type = self.tokens[idx].string
            possible_model_id = self.tokens[idx+2].string
            if case_fold_match:
                possible_model_id = possible_model_id.casefold()
            # the disambiguation model type cannot be Function
            if disambig_model_type == wc_lang.core.Function.__name__:
                return ("'{}', a {}.{}, contains '{}', which uses 'Function' as a disambiguation "
                        "model type but doesn't use Function syntax".format(self.expression, self.model_class.__name__,
                                                                            self.attribute, disambig_model_match))

            # the disambiguation model type must be in self.valid_models
            wc_lang_model_type = self.get_wc_lang_model_type(disambig_model_type)
            if wc_lang_model_type is None:
                return ("'{}', a {}.{}, contains '{}', but the disambiguation model type '{}' "
                        "cannot be referenced by '{}' expressions".format(
                            self.expression, self.model_class.__name__,
                            self.attribute, disambig_model_match, disambig_model_type,
                            self.model_class.__name__))

            if possible_model_id not in self._objects[wc_lang_model_type]:
                return "'{}', a {}.{}, contains '{}', but '{}' is not the id of a '{}'".format(
                    self.expression, self.model_class.__name__, self.attribute, disambig_model_match,
                    possible_model_id, disambig_model_type)

            return LexMatch([WcLangToken(TokCodes.wc_lang_obj_id, disambig_model_match, wc_lang_model_type,
                                         possible_model_id, self._objects[wc_lang_model_type][possible_model_id])],
                            len(self.model_type_disambig_pattern))

        # no match
        return None

    def related_object_id(self, idx, case_fold_match=False):
        """ Try to parse a related object `wc_lang` id from `self.tokens` at `idx`

        Different `wc_lang` objects match different Python token patterns. The default pattern
        is (token.NAME, ), but an object of type `model_type` can define a custom pattern in
        `model_type.Meta.token_pattern`, as Species does. Some patterns may consume multiple Python tokens.

        Args:
            idx (:obj:`int`): current index into `tokens`
            case_fold_match (:obj:`bool`, optional): if set, casefold identifiers before matching;
                identifier keys in `self._objects` must already be casefold'ed; default=False

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a related object id is found, return a `LexMatch` describing it.
        """
        token_matches = set()
        id_matches = set()
        for model_type in self.valid_models:
            token_pattern = (token.NAME, )
            if hasattr(model_type.Meta, 'token_pattern'):
                token_pattern = model_type.Meta.token_pattern
            match_string = self.match_tokens(token_pattern, idx)
            if match_string:
                token_matches.add(match_string)
                # is match_string the ID of an instance in model_type?
                if case_fold_match:
                    if match_string.casefold() in self._objects[model_type]:
                        id_matches.add(IdMatch(model_type, token_pattern, match_string))
                else:
                    if match_string in self._objects[model_type]:
                        id_matches.add(IdMatch(model_type, token_pattern, match_string))

        if not id_matches:
            if token_matches:
                return ("'{}', a {}.{}, contains the identifier(s) '{}', which aren't "
                        "the id(s) of an object".format(
                            self.expression, self.model_class.__name__,
                            self.attribute, "', '".join(token_matches)))
            return None

        if 1 < len(id_matches):
            # as lexers always do, pick the longest match
            id_matches_by_length = sorted(id_matches, key=lambda id_match: len(id_match.match_string))
            longest_length = len(id_matches_by_length[-1].match_string)
            longest_matches = set()
            while id_matches_by_length and len(id_matches_by_length[-1].match_string) == longest_length:
                longest_matches.add(id_matches_by_length.pop())
            id_matches = longest_matches

        if 1 < len(id_matches):
            # error: multiple, maximal length matches
            matches_error = ["'{}' as a {} id".format(id_val, model_type.__name__)
                             for model_type, _, id_val in sorted(id_matches, key=lambda id_match: id_match.model_type.__name__)]
            matches_error = ', '.join(matches_error)
            return "'{}', a {}.{}, contains multiple model object id matches: {}".format(
                self.expression, self.model_class.__name__, self.attribute, matches_error)

        else:
            # return a lexical match about a related id
            match = id_matches.pop()
            right_case_match_string = match.match_string
            if case_fold_match:
                right_case_match_string = match.match_string.casefold()
            return LexMatch(
                [WcLangToken(TokCodes.wc_lang_obj_id, match.match_string, match.model_type, right_case_match_string,
                             self._objects[match.model_type][right_case_match_string])],
                len(match.token_pattern))

    def fun_call_id(self, idx, case_fold_match='unused'):
        """ Try to parse a Python math function call from `self.tokens` at `idx`

        Each `wc_lang` object `model_class` that contains an expression which can use Python math
        functions must define the set of allowed functions in `Meta.valid_functions` of the
        model_class Expression Model.

        Args:
            idx (:obj:`int`): current index into `self.tokens`
            case_fold_match (:obj:`str`, optional): ignored keyword; makes `ParsedExpression.tokenize()` simpler

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a function call is found, return a `LexMatch` describing it.
        """
        fun_match = self.match_tokens(self.function_pattern, idx)
        if fun_match:
            fun_name = self.tokens[idx].string
            # function_pattern is "identifier ("
            # the closing paren ")" will simply be encoded as a WcLangToken with tok_code == op

            # are Python math functions defined?
            if not hasattr(self.model_class.Meta, 'valid_functions'):
                return ("'{}', a {}.{}, contains the func name '{}', but {}.Meta doesn't "
                        "define 'valid_functions'".format(self.expression,
                                                          self.model_class.__name__, self.attribute, fun_name, self.model_class.__name__))

            function_ids = set([f.__name__ for f in self.model_class.Meta.valid_functions])

            # is the function allowed?
            if fun_name not in function_ids:
                return ("'{}', a {}.{}, contains the func name '{}', but it isn't in "
                        "{}.Meta.valid_functions: {}".format(self.expression, self.model_class.__name__,
                                                             self.attribute, fun_name, self.model_class.__name__, ', '.join(function_ids)))

            # return a lexical match about a math function
            return LexMatch(
                [WcLangToken(TokCodes.math_fun_id, fun_name), WcLangToken(TokCodes.op, '(')],
                len(self.function_pattern))

        # no match
        return None

    def tokenize(self, case_fold_match=False):
        """ Tokenize a Python expression in `self.expression`

        Args:
            case_fold_match (:obj:`bool`, optional): if set, casefold identifiers before matching;
                identifier keys in `self._objects` must already be casefold'ed; default = False

        Returns:
            (:obj:`tuple`): either `(:obj:`None`, :obj:`list` of :obj:`str)` containing a list of errors, or
                `(:obj:`list`, :obj:`dict`)` containing a list of :obj:`WcLangToken`s and a dict of Model
                instances used by this list, grouped by Model type

        Raises:
            (:obj:`ParsedExpressionError`): if `model_class` does not have a `Meta` attribute
        """
        self.__reset_tokenization()

        # detect and report bad tokens
        bad_tokens = set()
        for tok in self.tokens:
            if tok.exact_type in self.illegal_tokens:
                if tok.string and tok.string != ' ':
                    bad_tokens.add(tok.string)
                else:
                    bad_tokens.add(token.tok_name[tok.type])
        if bad_tokens:
            self.errors.append("'{}', a {}.{}, contains bad token(s): '{}'".format(
                self.expression, self.model_class.__name__,
                self.attribute, "', '".join(bad_tokens)))
            return (None, self.errors)

        idx = 0
        while idx < len(self.tokens):

            # categorize token codes
            token_code = TokCodes.other
            if self.tokens[idx].type == token.OP:
                token_code = TokCodes.op
            elif self.tokens[idx].type == token.NUMBER:
                token_code = TokCodes.number

            # a token that isn't an identifier needs no processing
            if self.tokens[idx].type != token.NAME:
                # record non-identifier token
                self.wc_tokens.append(WcLangToken(token_code, self.tokens[idx].string))
                idx += 1
                continue

            matches = []
            tmp_errors = []
            for get_wc_lang_lexical_element in [self.related_object_id, self.disambiguated_id, self.fun_call_id]:
                result = get_wc_lang_lexical_element(idx, case_fold_match=case_fold_match)
                if result is not None:
                    if isinstance(result, str):
                        tmp_errors.append(result)
                    elif isinstance(result, LexMatch):
                        matches.append(result)
                    else:   # pragma no cover
                        assert False, "result is neither str nor LexMatch '{}'".format(result)

            # should find either matches or errors
            assert matches or tmp_errors, "No matches or errors found in '{}'".format(self.expression)
            # if only errors are found, break to return them
            if tmp_errors and not matches:
                self.errors = tmp_errors
                break

            # matches is a list of LexMatch, if it contains one longest match, use that, else report error
            # sort matches by Python token pattern length
            matches_by_length = sorted(matches, key=lambda lex_match: lex_match.num_py_tokens)
            longest_length = matches_by_length[-1].num_py_tokens
            longest_matches = []
            while matches_by_length and matches_by_length[-1].num_py_tokens == longest_length:
                longest_matches.append(matches_by_length.pop())
            assert len(longest_matches) <= 1, "multiple longest matches: '{}'".format(longest_matches)

            # good match
            # advance idx to the next token
            # record match data in self.wc_tokens and self.related_objects
            match = longest_matches.pop()
            idx += match.num_py_tokens
            wc_lang_tokens = match.wc_lang_tokens
            self.wc_tokens.extend(wc_lang_tokens)
            for wc_lang_token in wc_lang_tokens:
                if wc_lang_token.tok_code == TokCodes.wc_lang_obj_id:
                    self.related_objects[wc_lang_token.model_type][wc_lang_token.model_id] = wc_lang_token.model

        # detect ambiguous tokens
        valid_function_names = [func.__name__ for func in self.valid_functions]
        for wc_token in self.wc_tokens:
            if wc_token.tok_code in [TokCodes.wc_lang_obj_id, TokCodes.math_fun_id]:
                matching_items = []

                for model_type in self.valid_models:
                    if wc_token.token_string in self._objects[model_type]:
                        matching_items.append(model_type.__name__)

                if wc_token.token_string in valid_function_names:
                    matching_items.append('function')

                if len(matching_items) > 1:
                    self.errors.append('Token `{}` is ambiguous. Token matches a {} and a {}.'.format(
                        wc_token.token_string, ', a '.join(matching_items[0:-1]), matching_items[-1]))

        if self.errors:
            return (None, self.errors)
        return (self.wc_tokens, self.related_objects)

    def test_eval(self, test_val=1.0):
        """ Test evaluate the expression with the value of all used models equal to `test_val`.

        Called to validate this `ParsedExpression`.

        Args:
            test_val (:obj:`float`, optional): the value assumed for used Models

        Returns:
            :obj:`object`: the value of the expression

        Raises:
            :obj:`ParsedExpressionError`: if the expression evaluation fails
        """
        # do not eval an expression that could not be tokenized
        if not self.wc_tokens:
            raise ParsedExpressionError("cannot evaluate '{}', as it not been successfully tokenized".format(
                self.expression))

        model_vals = {}
        for model_type, models in self.related_objects.items():
            model_vals[model_type] = {}
            for model in models.values():
                model_vals[model_type][model] = test_val

        return self.eval(model_vals)

    def eval(self, model_vals):
        """ Evaluate the expression

        This expression must have been successfully `tokenize`d.

        Approach:
            * Replace references to used Models in `self.wc_tokens` with `test_val`
            * Join the elements of `self.wc_tokens` into a Python expression
            * `eval` the Python expression

        Args:
            model_vals (:obj:`dict` of `Model): dictionary of values of the models; maps
                model type to a dict mapping models to Model values

        Returns:
            :obj:`object`: the value of the expression

        Raises:
            :obj:`ParsedExpressionError`: if the expression evaluation fails
        """
        # do not eval an expression that could not be tokenized
        if not self.wc_tokens:
            raise ParsedExpressionError("cannot evaluate '{}', as it not been successfully tokenized".format(
                self.expression))

        evaled_tokens = []
        idx = 0
        while idx < len(self.wc_tokens):
            wc_token = self.wc_tokens[idx]
            if wc_token.tok_code == TokCodes.wc_lang_obj_id:
                evaled_tokens.append(str(model_vals[wc_token.model_type][wc_token.model]))
                if wc_token.model_type == wc_lang.core.Function:
                    # skip past the following ( ) tokens -- they're just syntactic sugar for Functions
                    idx += 2
            else:
                evaled_tokens.append(wc_token.token_string)
            idx += 1

        expression = ' '.join(evaled_tokens)
        local_ns = {func.__name__: func for func in self.valid_functions}

        error_suffix = " cannot eval expression '{}' in {}; ".format(expression,
                                                                     self.model_class.__name__)

        try:
            return eval(expression, {}, local_ns)
        except SyntaxError as error:
            raise ParsedExpressionError("SyntaxError:" + error_suffix + str(error))
        except NameError as error:  # pragma: no cover
            raise ParsedExpressionError("NameError:" + error_suffix + str(error))
        except Exception as error:  # pragma: no cover
            raise ParsedExpressionError("Exception:" + error_suffix + str(error))

    def __str__(self):
        rv = []
        rv.append("model_class: {}".format(self.model_class.__name__))
        rv.append("expression: '{}'".format(self.expression))
        rv.append("attribute: {}".format(self.attribute))
        rv.append("tokens: {}".format("'"+"', '".join([t.string for t in self.tokens])+"'"))
        rv.append("related_objects: {}".format(self.related_objects))
        rv.append("errors: {}".format(self.errors))
        rv.append("wc_tokens: {}".format(self.wc_tokens))
        return '\n'.join(rv)


class ExpressionVerifier(object):
    """ Verify whether a sequence of `WcLangToken` tokens

    An `ExpressionVerifier` consists of two parts:

    * An optional method `valid_wc_lang_tokens` that examines the content of individual tokens
      and returns `(True, True)` if they are all valid, or (`False`, error) otherwise. It can be
      overridden by subclasses.
    * A `DFSMAcceptor` which determines whether the tokens describe a particular pattern
      `validate()` combines these parts.

    Attributes:
        dfsm_acceptor (:obj:`DFSMAcceptor`): the DFSM acceptor
        empty_is_valid (:obj:`bool`): if set, then an empty sequence of tokens is valid
    """

    def __init__(self, start_state, accepting_state, transitions, empty_is_valid=False):
        self.dfsm_acceptor = DFSMAcceptor(start_state, accepting_state, transitions)
        self.empty_is_valid = empty_is_valid

    def valid_wc_lang_tokens(self, wc_lang_tokens):
        return (True, True)

    def make_dfsa_messages(self, wc_lang_tokens):
        """ Convert a sequence of `WcLangToken`s into a list of messages for transitions

        Args:
            wc_lang_tokens (:obj:`iterator` of `WcLangToken`): sequence of `WcLangToken`s

        Returns:
            :obj:`object`: `None` if `wc_lang_tokens` cannot be converted into a sequence of messages,
                or a `list` of `tuple` of pairs (token code, message modifier)
        """
        messages = []
        for wc_lang_tok in wc_lang_tokens:
            messages.append((wc_lang_tok.tok_code, None))
        return messages

    def validate(self, wc_lang_tokens):
        """ Indicate whether `wc_lang_tokens` is valid

        Args:
            wc_lang_tokens (:obj:`iterator` of `WcLangToken`): sequence of `WcLangToken`s

        Returns:
            :obj:`tuple`: (`False`, error) if `wc_lang_tokens` is valid, or (`True`, `None`) if it is
        """
        if self.empty_is_valid and not wc_lang_tokens:
            return (True, None)
        valid, error = self.valid_wc_lang_tokens(wc_lang_tokens)
        if not valid:
            return (False, error)
        dfsa_messages = self.make_dfsa_messages(wc_lang_tokens)
        if DFSMAcceptor.ACCEPT == self.dfsm_acceptor.run(dfsa_messages):
            return (True, None)
        else:
            return (False, "Not a linear expression")


class LinearExpressionVerifier(ExpressionVerifier):
    """ Verify whether a sequence of tokens (`WcLangToken`s) describes a linear function of identifiers

    In particular, a valid linear expression must have the structure:
        * `(identifier | number '*' identifier) (('+' | '-') (identifier | number '*' identifier))*`
    """

    # Transitions in valid linear expression
    TRANSITIONS = [   # (current state, message, next state)
        ('need number or id', (TokCodes.number, None), 'need * id'),
        ('need * id', (TokCodes.op, '*'), 'need id'),
        ('need id', (TokCodes.wc_lang_obj_id, None), 'need + | - | end'),
        ('need number or id', (TokCodes.wc_lang_obj_id, None), 'need + | - | end'),
        ('need + | - | end', (TokCodes.op, '+'), 'need number or id'),
        ('need + | - | end', (TokCodes.op, '-'), 'need number or id'),
        ('need + | - | end', (None, None), 'end'),
    ]

    def __init__(self):
        super().__init__(start_state='need number or id', accepting_state='end',
                         transitions=self.TRANSITIONS, empty_is_valid=True)

    def valid_wc_lang_tokens(self, wc_lang_tokens):
        """ Check whether the content of a sequence of `WcLangToken`s is valid

        In particular, all numbers in `wc_lang_tokens` must be floats, and all token codes must not
        be `math_fun_id` or `other`.

        Args:
            wc_lang_tokens (:obj:`iterator` of `WcLangToken`): sequence of `WcLangToken`s

        Returns:
            :obj:`tuple`: (`False`, error) if `wc_lang_tokens` cannot be a linear expression, or
                (`True`, `True`) if it can
        """
        for wc_lang_tok in wc_lang_tokens:
            if wc_lang_tok.tok_code in set([TokCodes.math_fun_id, TokCodes.other]):
                return (False, "messages do not use token codes `math_fun_id` or `other`")
            if wc_lang_tok.tok_code == TokCodes.number:
                try:
                    float(wc_lang_tok.token_string)
                except ValueError as e:
                    return (False, str(e))

        return (True, True)

    def make_dfsa_messages(self, wc_lang_tokens):
        """ Convert a sequence of `WcLangToken`s into a list of messages for transitions in 
        `LinearExpressionVerifier.TRANSITIONS`

        Args:
            wc_lang_tokens (:obj:`iterator` of `WcLangToken`): sequence of `WcLangToken`s

        Returns:
            :obj:`object`: `None` if `wc_lang_tokens` cannot be converted into a sequence of messages
                to validate a linear expression, or a `list` of `tuple` of pairs (token code, message modifier)
        """
        messages = []
        for wc_lang_tok in wc_lang_tokens:
            message_tok_code = wc_lang_tok.tok_code
            if wc_lang_tok.tok_code == TokCodes.wc_lang_obj_id:
                message_modifier = None
            elif wc_lang_tok.tok_code == TokCodes.number:
                message_modifier = None
            elif wc_lang_tok.tok_code == TokCodes.op:
                message_modifier = wc_lang_tok.token_string
            else:
                return None
            messages.append((message_tok_code, message_modifier))
        messages.append((None, None))
        return messages
