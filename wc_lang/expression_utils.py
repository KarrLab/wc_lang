'''Static methods which handle rate laws and their expressions

:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2017-10-23
:Copyright: 2016-2017, Karr Lab
:License: MIT
'''
import tokenize
import token
from six.moves import cStringIO
from io import BytesIO
from math import isnan, ceil, floor, exp, pow, log, log10
from collections import namedtuple

from wc_utils.util.list import det_dedupe
from wc_utils.util.enumerate import CaseInsensitiveEnum
import wc_lang


CONCENTRATIONS_DICT = 'concentrations'
PARAMETERS_DICT = 'parameters'


# TODOS
'''
build
    support ModelName. disambiguation prefixes to IDs
    support wc_lang.Function instances
    make dynamic_model.eval_dynamic_obj() handle all dyamic Models
    create a WcLangExpressionError, and replace ValueError with it
    make generic ModelWithExpression class in wc_lang
    add test_eval() method
    test with real WC models
    Jupyter example
    replace existing RE parsing and expression eval code
    use ExpressionUtils to deserialize all wc_lang expressions
    rename Function to Macro
cleanup
    have valid_functions defined as sets, not tuples
    stop using & and remove RateLawUtils
'''
class RateLawUtils(object):
    '''A set of static rate_law methods '''

    @staticmethod
    def transcode(rate_law_equation, species, parameters):
        '''Translate a `wc_lang.core.RateLawEquation` into a python expression that can be evaluated
        during a simulation.

        Args:
            rate_law_equation (:obj:`wc_lang.core.RateLawEquation`): a rate law equation
            species (:obj:`set` of `wc_lang.core.Species`): the species that use the rate law
            parameters (:obj:`set` of `wc_lang.core.Parameter`): the parameters that use the rate law

        Returns:
            The python expression, or None if the rate law doesn't have an equation

        Raises:
            ValueError: If `rate_law_equation` contains `__`, which increases its security risk, or
            if `rate_law_equation` refers to species not in `species` or parameters not in `parameters`
        '''
        def possible_specie_id(tokens):
            '''Indicate whether `tokens` begins with 4 tokens whose syntax id_matches a specie ID

            Specie IDs have the form `string[string]`, as documented in `wc_lang.core.Species.id()`.

            Args:
                tokens (:obj:`list` of (token_num, token_val)): a list of Python tokens

            Returns:
                :obj:`bool`: True if the initial elements of `tokens` has the syntax of a specie ID
            '''
            if len(tokens) < 4:
                return False
            toknums = [token_tmp[0] for token_tmp in tokens]
            specie_id_pattern = [token.NAME, token.OP, token.NAME, token.OP]
            if toknums == specie_id_pattern:
                tokvals = [token_tmp[1] for token_tmp in tokens]
                return tokvals[1] == '[' and tokvals[3] == ']'
            return False

        def possible_parameter_id(tokens, reserved_tokens):
            '''Indicate whether `tokens` begins with 4 tokens whose syntax id_matches a parameter ID

            Parameter IDs have the form `string`, as documented in `wc_lang.core.Parameter.id`.

            Args:
                tokens (:obj:`list` of (token_num, token_val)): a list of Python tokens

            Returns:
                :obj:`bool`: True if the initial elements of `tokens` has the syntax of a parameter ID
            '''
            return len(tokens) >= 1 and \
                tokens[0][0] == token.NAME and \
                tokens[0][1] not in reserved_tokens

        def convert_specie_name(tokens, species_ids, rate_law_expression):
            '''Translate a specie ID in `tokens` into a Python expression to be eval'ed during simulation

            Args:
                tokens (:obj:`list` of (token_num, token_val)): a list of 4 Python tokens that
                    have the syntax of a specie ID
                species_ids (:obj:`set`): IDs of the species used by the rate law expression
                rate_law_expression (:obj:`string`): the rate law expression being transcoded

            Returns:
                (:obj:`string`): a Python expression, transcoded to look up the specie concentration
                    in `concentrations[]`

            Raises:
                ValueError: if `tokens` does not represent a specie in `species_ids`
            '''
            tokvals = [token_tmp[1] for token_tmp in tokens]
            parsed_id = wc_lang.Species.gen_id(tokvals[0], tokvals[2])
            if parsed_id in species_ids:
                return " {}['{}']".format(CONCENTRATIONS_DICT, parsed_id)
            else:
                raise ValueError("'{}' not a known specie in rate law '{}'".format(
                    parsed_id, rate_law_expression))

        def convert_parameter_name(tokens, parameter_ids, rate_law_expression):
            '''Translate a parameter ID in `tokens` into a Python expression to be eval'ed during simulation

            Args:
                tokens (:obj:`list` of (token_num, token_val)): a list of 1 Python token that
                    have the syntax of a parameter ID
                parameter_ids (:obj:`set`): IDs of the parameter used by the rate law expression
                rate_law_expression (:obj:`string`): the rate law expression being transcoded

            Returns:
                (:obj:`string`): a Python expression, transcoded to look up the parameter value
                    in `parameters[]`

            Raises:
                ValueError: if `tokens` does not represent a parameter in `parameter_ids`
            '''
            parsed_id = tokens[0][1]
            if parsed_id in parameter_ids:
                return " {}['{}']".format(PARAMETERS_DICT, parsed_id)
            else:
                raise ValueError("'{}' not a known parameter in rate law '{}'".format(
                    parsed_id, rate_law_expression))

        if '__' in rate_law_equation.expression:
            raise ValueError("Security risk: rate law expression '{}' contains '__'.".format(
                rate_law_equation.expression))

        rate_law_expression = rate_law_equation.expression
        species_ids = set([specie.id() for specie in species])
        parameter_ids = set([parameter.id for parameter in parameters])
        reserved_parameter_ids = set(map(lambda f: f.__name__, wc_lang.RateLawEquation.Meta.valid_functions)) | \
            set(['k_cat', 'k_m'])

        # Rate laws must be tokenized to properly construct a Python expression.
        # A prior implementation which used REs and string replace() contained a subtle bug that
        # was triggered when one species name matched the suffix of another.
        # For example consider a rate law with
        #   expression = 'xy[c] + y[c]'
        # Replacing the species 'y[c]' in this expression with "concentration['y[c]']" like this
        #   expression.replace('y[c]', "concentration['y[c]']")
        # would erroneously generate
        #   "xconcentration['y[c]'] + concentration['y[c]']"
        g = tokenize.generate_tokens(cStringIO(rate_law_equation.expression).readline)
        tokens = [(toknum, tokval) for toknum, tokval, _, _, _ in g]
        result = []
        idx = 0
        while idx < len(tokens):
            if possible_specie_id(tokens[idx:idx+4]):
                result.append(
                    (token.NAME, convert_specie_name(tokens[idx:idx+4], species_ids, rate_law_expression)))
                idx += 4
            elif possible_parameter_id(tokens[idx:idx+1], reserved_parameter_ids):
                result.append(
                    (token.NAME, convert_parameter_name(tokens[idx:idx+1], parameter_ids, rate_law_expression)))
                idx += 1
            else:
                result.append((tokens[idx]))
                idx += 1
        py_expression = tokenize.untokenize(result)
        return py_expression

    @staticmethod
    def transcode_rate_laws(model):
        '''Transcode all of `model`'s rate law expressions into Python expressions

        Args:
            model (:obj:`wc_lang.core.Model`): a `wc_lang.core.Model`

        Raises:
            ValueError: If a rate law cannot be transcoded
        '''
        for submodel in model.get_submodels():
            for reaction in submodel.reactions:
                for rate_law in reaction.rate_laws:
                    if hasattr(rate_law, 'equation'):
                        rate_law.equation.transcoded = RateLawUtils.transcode(rate_law.equation,
                                                                              submodel.get_species(),
                                                                              rate_law.equation.parameters)

    @staticmethod
    def eval_reaction_rate_laws(reaction, concentrations, parameters):
        '''Evaluate a reaction's rate laws at the given species concentrations and parameter values

        Args:
            reaction (:obj:`wc_lang.core.Reaction`): a Reaction instance
            concentrations (:obj:`dict` of :obj:`species_id` -> `float`):
                a dictionary of species concentrations
            parameters (:obj:`dict` of :obj:`parameter_id` -> `float`):
                a dictionary of parameter values

        Returns:
            (:obj:`list` of `float`): the reaction's forward and, perhaps, backward rates
        '''
        rates = []
        for rate_law in reaction.rate_laws:
            rates.append(RateLawUtils.eval_rate_law(rate_law, concentrations, parameters))
        return rates

    @staticmethod
    def eval_rate_law(rate_law, concentrations, parameters, transcoded_equation=None):
        '''Evaluate a rate law at the given species concentrations and parameter values

        Args:
            rate_law (:obj:`wc_lang.core.RateLaw`): a RateLaw instance
            concentrations (:obj:`dict` of :obj:`species_id` -> `float`):
                a dictionary of species concentrations
            parameters (:obj:`dict` of :obj:`parameter_id` -> `float`):
                a dictionary of parameter values
            transcoded_equation (:obj:`str`, optional): the rate law's transcoded_equation; if not
                provided, will be taken from rate_law.equation.transcoded

        Returns:
            :obj:`float`): the rate law's rate

        Raises:
            ValueError: if the rate law has a syntax error
            NameError: if the rate law references a specie whose concentration is not provided in
                `concentrations` or a parameter whose value is not provided in `parameters`
            Exception: if the rate law has other errors, such as a reference to an unknown function
        '''
        if not transcoded_equation:
            transcoded_equation = rate_law.equation.transcoded

        local_ns = {func.__name__: func for func in wc_lang.RateLawEquation.Meta.valid_functions}
        if hasattr(rate_law, 'k_cat') and not isnan(rate_law.k_cat):
            local_ns['k_cat'] = rate_law.k_cat
        if hasattr(rate_law, 'k_m') and not isnan(rate_law.k_m):
            local_ns['k_m'] = rate_law.k_m
        local_ns[CONCENTRATIONS_DICT] = concentrations
        local_ns[PARAMETERS_DICT] = parameters

        try:
            return eval(transcoded_equation, {}, local_ns)
        except SyntaxError as error:
            raise ValueError("Error: syntax error in transcoded rate law '{}'.".format(transcoded_equation))
        except NameError as error:
            raise NameError("Error: NameError in transcoded rate law '{}': '{}'".format(
                transcoded_equation, error))
        except Exception as error:
            raise Exception("Error: unable to eval transcoded rate law '{}': {} for {}".format(
                transcoded_equation, error.__class__.__name__, error))


class TokCodes(int, CaseInsensitiveEnum):
    """ Token codes used in parsed expressions """
    wc_lang_obj_id = 1
    math_fun_id = 2
    other = 3


# result returned by a tokens lexer, like disambiguated_id()
LexMatch = namedtuple('LexMatch', 'wc_lang_tokens, num_py_tokens')
LexMatch.__doc__ += ': result returned by a lexer method that matches a wc_lang expression element'
LexMatch.wc_lang_tokens.__doc__ = 'List of WcLangTokens created'
LexMatch.num_py_tokens.__doc__ = 'Number of Python tokens consumed'


# a matched token pattern used by deserialize
IdMatch = namedtuple('IdMatch', 'model_type, token_pattern, match_string')
IdMatch.__doc__ += ': Matched token pattern used by deserialize'
IdMatch.model_type.__doc__ = 'The type of Model matched'
IdMatch.token_pattern.__doc__ = 'The token pattern used by the match'
IdMatch.match_string.__doc__ = 'The matched string'


# a token in a parsed wc_lang expression, returned in a list by deserialize
WcLangToken = namedtuple('WcLangToken', 'tok_code, token_string, model_type, model_id')
# make model_type and model_id optional: see https://stackoverflow.com/a/18348004
WcLangToken.__new__.__defaults__ = (None, None)
WcLangToken.__doc__ += ': Token in a parsed wc_lang expression'
WcLangToken.tok_code.__doc__ = 'TokCodes encoding'
WcLangToken.token_string.__doc__ = "The token's string"
WcLangToken.model_type.__doc__ = "When tok_code is wc_lang_obj_id, the wc_lang obj's type"
WcLangToken.model_id.__doc__ = "When tok_code is wc_lang_obj_id, the wc_lang obj's id"


class WcLangExpression(object):
    """ An expression in a wc_lang Model

    Attributes:
        model_class (:obj:`obj_model.Model`): the `wc_lang` `Model` which has an expression
        attribute (:obj:`str`): the attribute name of the expression in `model_class`
        expression (:obj:`str`): the expression defined in the wc_lang Model
        tokens (:obj:`list` of :obj:`namedtuple`): a list of tokens generated by `tokenize.tokenize()`
        idx (:obj:`int`): current index into `tokens`
        objects (:obj:`dict`): dict of wc_lang Models that might be referenced in expression; maps
            model type to set of ids for the object type
        related_objects (:obj:`dict`): models that are referenced in `expression`; maps model type to
            list of ids used
        errors (:obj:`list` of `str`): errors found when parsing an `expression` fails
        wc_tokens (:obj:`list` of `WcLangToken`): tokens obtained when an `expression` is successfully parsed
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
        """ Create an instance of WcLangExpression

        Raises:
            (:obj:`ValueError`): if lexical analysis of `expression` raises an exception
        """
        self.model_class = model_class
        self.attribute = attribute
        # strip leading and trailing whitespace from expression, which would create a bad token error
        self.expression = expression.strip()
        self.objects = objects

        self.related_objects = {}
        for model_type in objects.keys():
            self.related_objects[model_type] = []

        try:
            g = tokenize.tokenize(BytesIO(self.expression.encode('utf-8')).readline)
            # strip the leading ENCODING and trailing ENDMARKER tokens
            self.tokens = list(g)[1:-1]
        except tokenize.TokenError as e:
            raise ValueError("parsing '{}', a {}.{}, creates a Python syntax error: '{}'".format(
                self.expression, self.model_class.__name__, self.attribute, str(e)))

        self.idx = 0
        self.errors = []
        self.wc_tokens = []

    def get_wc_lang_model_type(self, model_type_name):
        """ Find the `wc_lang` model type corresponding to `model_type_name`

        Args:
            model_type_name (:obj:`str`): the name of a purported `wc_lang` model type in an expression

            Returns:
                :obj:`object`: `None` if no model named `model_type_name` exists in `self.objects`,
                    else the type of the model with that name
        """
        for model_type in self.objects.keys():
            if model_type_name == model_type.__name__:
                return model_type
        return None

    def match_tokens(self, token_pattern, idx):
        """ Indicate whether `tokens` begins with a pattern of tokens that match `token_pattern`

        Args:
            token_pattern (:obj:`tuple` of `int`): a tuple of Python token numbers, taken from the
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
            if 0<tok_idx and self.tokens[idx+tok_idx-1].end != self.tokens[idx+tok_idx].start:
                return False
        match_val = ''.join([self.tokens[idx+i].string for i in range(len(token_pattern))])
        return match_val

    def disambiguated_id(self, idx):
        """ Try to parse a disambuated `wc_lang` id from `self.tokens` at `idx`

        Look for a disambugated id (either a Function written as `Function.identifier()`, or a
        Model written as `ModelType.model_id`). If tokens do not match, return `None`. If tokens match,
        but their values are wrong, return an error `str`.
        If a disambugated id is found, return a `LexMatch` describing it.

        Args:
            idx (:obj:`int`): current index into `tokens`

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a disambugated id is found, return a `LexMatch` describing it.
        """
        fun_match = self.match_tokens(self.fun_type_disambig_patttern, idx)
        if fun_match:
            purported_macro_name = self.tokens[idx+2].string
            # the disambiguation model type must be Function
            if self.tokens[idx].string != wc_lang.core.Function.__name__:
                return ("'{}', a {}.{}, contains '{}', which doesn't use 'Function' as a disambiguation "
                    "model type".format(self.expression, self.model_class.__name__, self.attribute, fun_match))
            # the identifier must be in the Function objects
            if wc_lang.core.Function not in self.objects or purported_macro_name not in self.objects[wc_lang.core.Function]:
                return "'{}', a {}.{}, contains '{}', which doesn't refer to a Function in 'objects'".format(
                    self.expression, self.model_class.__name__, self.attribute, fun_match)
            return LexMatch([WcLangToken(TokCodes.wc_lang_obj_id, fun_match, wc_lang.core.Function, purported_macro_name)],
                len(self.fun_type_disambig_patttern))

        disambig_model_match = self.match_tokens(self.model_type_disambig_pattern, idx)
        if disambig_model_match:
            disambig_model_type = self.tokens[idx].string
            purported_model_name = self.tokens[idx+2].string
            # the disambiguation model type cannot be Function
            if disambig_model_type == wc_lang.core.Function.__name__:
                return ("'{}', a {}.{}, contains '{}', which uses 'Function' as a disambiguation "
                    "model type but doesn't use Function syntax".format(self.expression, self.model_class.__name__,
                    self.attribute, disambig_model_match))

            # the disambiguation model type must be in self.objects
            wc_lang_model_type = self.get_wc_lang_model_type(disambig_model_type)
            if wc_lang_model_type is None:
                return ("'{}', a {}.{}, contains '{}', but the disambiguation model type '{}' "
                    "cannot be referenced by '{}' expressions".format(self.expression, self.model_class.__name__,
                    self.attribute, disambig_model_match, disambig_model_type, self.model_class.__name__))

            if purported_model_name not in self.objects[wc_lang_model_type]:
                return "'{}', a {}.{}, contains '{}', but '{}' is not the id of a '{}'".format(
                    self.expression, self.model_class.__name__, self.attribute, disambig_model_match,
                    purported_model_name, disambig_model_type)

            return LexMatch([WcLangToken(TokCodes.wc_lang_obj_id, disambig_model_match, wc_lang_model_type, purported_model_name)],
                len(self.model_type_disambig_pattern))

        # no match
        return None

    def related_object_id(self, idx):
        """ Try to parse a related object `wc_lang` id from `self.tokens` at `idx`

        Different `wc_lang` objects match different Python token patterns. The default pattern
        is (token.NAME, ), but an object of type `model_type` can define a custom pattern in
        `model_type.Meta.token_pattern`, as Species does. Some patterns may consume multiple Python tokens.

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a related object id is found, return a `LexMatch` describing it.
        """
        token_matches = set()
        id_matches = set()
        object_type_names = [model_type.__name__ for model_type in self.objects.keys()]
        for model_type in self.objects.keys():
            token_pattern = (token.NAME, )
            if hasattr(model_type.Meta, 'token_pattern'):
                token_pattern = model_type.Meta.token_pattern
            match_string = self.match_tokens(token_pattern, idx)
            if match_string:
                token_matches.add(match_string)
                # is match_string the ID of an instance in model_type?
                if match_string in self.objects[model_type]:
                    id_matches.add(IdMatch(model_type, token_pattern, match_string))

        if not id_matches:
            if token_matches:
                return ("'{}', a {}.{}, contains the identifier(s) '{}', which aren't "
                    "the id(s) of an object in 'objects'".format(self.expression,
                    self.model_class.__name__, self.attribute, "', '".join(token_matches)))
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
            # error: multiple, maximal length model id_matches
            matches_error = ["'{}' as a {} id".format(id_val, model_type.__name__)
                for model_type, _, id_val in sorted(id_matches, key=lambda id_match: id_match.model_type.__name__)]
            matches_error = ', '.join(matches_error)
            return "'{}', a {}.{}, contains multiple model object id_matches: {}".format(self.expression,
                self.model_class.__name__, self.attribute, matches_error)

        else:
            # return a lexical match about a related id
            match = id_matches.pop()
            return LexMatch(
                [WcLangToken(TokCodes.wc_lang_obj_id, match.match_string, match.model_type, match.match_string)],
                len(match.token_pattern))

    def fun_call_id(self, idx):
        """ Try to parse a Python math function call from `self.tokens` at `idx`

        Each `wc_lang` object `model_class` that contains an expression which can use Python math
        functions must define the set of allowed functions in `model_class.Meta.valid_functions`.

        Args:
            idx (:obj:`int`): current index into `self.tokens`

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a function call is found, return a `LexMatch` describing it.
        """
        fun_match = self.match_tokens(self.function_pattern, idx)
        if fun_match:
            fun_name = self.tokens[idx].string
            # function_pattern is "identifier ("
            # the closing paren ")" will simply be encoded as a WcLangToken with tok_code == other

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
                [WcLangToken(TokCodes.math_fun_id, fun_name), WcLangToken(TokCodes.other, '(')],
                len(self.function_pattern))

        # no match
        return None

    def lex_input(idx):
        matches = []
        errors = []
        for get_tok_fun in []:
            result = get_tok_fun(idx)
            if result is not None:
                if isinstance(result, str):
                    errors.append(result)
                elif isinstance(result, LexMatch):
                    matches.append(result)
        if not matches:
            pass
            # return all errors
        # sort matches by length
        # if one longest match, use it
        # else error

    def __str__(self):
        rv = []
        rv.append("model_class: {}".format(self.model_class.__name__))
        rv.append("expression: '{}'".format(self.expression))
        rv.append("attribute: {}".format(self.attribute))
        rv.append("tokens: {}".format("'"+"', '".join([t.string for t in self.tokens])+"'"))
        rv.append("idx: {}".format(self.idx))
        rv.append("objects: {}".format(self.objects))
        rv.append("related_objects: {}".format(self.related_objects))
        rv.append("errors: {}".format(self.errors))
        rv.append("wc_tokens: {}".format(self.wc_tokens))
        return '\n'.join(rv)


class ExpressionUtils(object):
    """ Utilities for processing expressions in Models

    Expressions are currently (July, 2018) used in four `wc_lang` `Model`s: `RateLawEquation`, `Function`,
    `StopCondition` (which is just a special case of `Function` that returns a boolean), and `ObjectiveFunction`.
    These expressions are limited Python expressions with specific semantics:

    * They must be syntactically correct Python.
    * No Python keywords, strings, or tokens that do not belong in expressions are allowed.
    * All Python identifiers must be the ID of an object in a whole-cell model, or components of
        an ID of an object in the model, or the name of a function in the `math` package. Objects in the model
        are provided in `objects`, and the allowed subset of functions in `math` must be provided in an
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
    """

    function_pattern = (token.NAME, token.LPAR)

    # enumerate and detect Python tokens that are illegal in wc_lang expressions
    # TODO: consider a handful of other tokens may also be illegal: COMMA, SEMI, TILDE, CIRCUMFLEX, and AT
    illegal_token_names = ['ENDMARKER', 'NEWLINE', 'INDENT', 'DEDENT', 'COLON', 'LBRACE', 'RBRACE',
        'PLUSEQUAL', 'MINEQUAL', 'STAREQUAL', 'SLASHEQUAL', 'PERCENTEQUAL', 'AMPEREQUAL', 'VBAREQUAL',
        'CIRCUMFLEXEQUAL', 'LEFTSHIFTEQUAL', 'RIGHTSHIFTEQUAL', 'DOUBLESTAREQUAL', 'DOUBLESLASHEQUAL',
        'ATEQUAL', 'RARROW', 'ELLIPSIS', 'AWAIT', 'ASYNC', 'ERRORTOKEN', 'N_TOKENS', 'NT_OFFSET']
    illegal_tokens = set()
    for illegal_name in illegal_token_names:
        illegal_tokens.add(getattr(token, illegal_name))

    @staticmethod
    def match_tokens(token_pattern, tokens):
        """ Indicate whether `tokens` begins with a pattern of tokens that match `token_pattern`

        Args:
            token_pattern (:obj:`tuple` of `int`): a tuple of Python token numbers, taken from the
            `token` module
            tokens (:obj:`list` of :obj:`namedtuple`): a list of tokens generated by
            `tokenize.tokenize()`

            Returns:
                :obj:`object`: :obj:`bool`, False if the initial elements of `tokens` do not match the
                syntax in `token_pattern`, or :obj:`str`, the matching string
        """
        if not token_pattern:
            return False
        if len(tokens) < len(token_pattern):
            return False
        for idx, token_pat_num in enumerate(token_pattern):
            if tokens[idx].exact_type != token_pat_num:
                return False
            # because a wc_lang ID shouldn't contain white space, do not allow it between the tokens
            # that match token_pattern
            if 0<idx and tokens[idx-1].end != tokens[idx].start:
                return False
        match_val = ''.join([tokens[i].string for i in range(len(token_pattern))])
        return match_val

    @staticmethod
    def deserialize(model_class, attribute, expression, objects):
        """ Deserialize a Python expression in `expression`

        Args:
            model_class (:obj:`obj_model.Model`): the `wc_lang` `Model` whose expression is being deserialize
            attribute (:obj:`str`): the attribute name of the expression in `model_class`
            expression (:obj:`str`): a `wc_lang` Python expression
            objects (:obj:`dict` of `list` of `obj_model.Model`): map from `obj_model.Model`s to
                `list`s of their instances

        Returns:
            (:obj:`tuple`): either `(None, :obj:list of str)` containing a list of errors, or `(:obj:list, :obj:dict)`
            containing a list of `WcLangToken`s and a dict of modifiers used in this list

        Raises:
            (:obj:`ValueError`): if `model_class` does not have a `Meta` attribute
        """
        errors = []
        wc_tokens = []
        related_objects = {}
        for model_type in objects.keys():
            related_objects[model_type] = []

        # strip leading and trailing whitespace from expression, which would create a bad token error
        expression = expression.strip()

        g = tokenize.tokenize(BytesIO(expression.encode('utf-8')).readline)
        # strip the leading ENCODING and trailing ENDMARKER tokens
        tokens = list(g)[1:-1]

        # detect and report bad tokens
        bad_tokens = set()
        for tok in tokens:
            if tok.exact_type in ExpressionUtils.illegal_tokens:
                if tok.string and tok.string != ' ':
                    bad_tokens.add(tok.string)
                else:
                    bad_tokens.add(token.tok_name[tok.type])
        if bad_tokens:
            error = "'{}', a {}.{}, contains bad token(s): '{}'".format(expression, model_class.__name__,
                attribute, "', '".join(bad_tokens))
            return (None, [error])

        idx = 0
        while idx < len(tokens):

            # a token that isn't an identifier needs no processing
            if tokens[idx].type != token.NAME:
                # record non-identifier token
                wc_tokens.append(WcLangToken(TokCodes.other, tokens[idx].string))
                idx += 1
                continue

            # look for a function call
            fun_match = ExpressionUtils.match_tokens(ExpressionUtils.function_pattern, tokens[idx:])
            if fun_match:
                fun_name = tokens[idx].string
                # function_pattern is "identifier (" which can reference either a wc_lang Macro or a Python math function
                try:

                    # are Python math functions defined?
                    if not hasattr(model_class.Meta, 'valid_functions'):
                        errors.append("'{}', a {}.{}, contains the func name '{}', but {}.Meta doesn't "
                            "define 'valid_functions'".format(expression,
                            model_class.__name__, attribute, fun_name, model_class.__name__))
                        idx += len(ExpressionUtils.function_pattern)
                        continue
                    function_ids = set([f.__name__ for f in model_class.Meta.valid_functions])

                    # is the function allowed?
                    if fun_name not in function_ids:
                        errors.append("'{}', a {}.{}, contains the func name '{}', but it isn't in "
                            "{}.Meta.valid_functions: {}".format(expression, model_class.__name__,
                            attribute, fun_name, model_class.__name__, ', '.join(function_ids)))
                        idx += len(ExpressionUtils.function_pattern)
                        continue

                    # record function tokens
                    wc_tokens.append(WcLangToken(TokCodes.math_fun_id, fun_name))
                    # add the '(' which was used to find the function
                    wc_tokens.append(WcLangToken(TokCodes.other, '('))
                    idx += len(ExpressionUtils.function_pattern)
                    continue
                except AttributeError as e:
                    raise ValueError(str(e))

            # look for a related object id; some ids consume multiple Python tokens
            token_matches = set()
            id_matches = set()
            object_type_names = [model_type.__name__ for model_type in objects.keys()]
            for model_type in objects.keys():
                token_pattern = (token.NAME, )
                if hasattr(model_type.Meta, 'token_pattern'):
                    token_pattern = model_type.Meta.token_pattern
                match_string = ExpressionUtils.match_tokens(token_pattern, tokens[idx:])
                if match_string:
                    token_matches.add(match_string)
                    # is match_string the ID of an instance in model_type?
                    if match_string in objects[model_type]:
                        id_matches.add(IdMatch(model_type, token_pattern, match_string))

            if not id_matches:
                if token_matches:
                    errors.append("'{}', a {}.{}, contains the identifier(s) '{}', which aren't "
                        "the id(s) of an object in 'objects'".format(expression,
                        model_class.__name__, attribute, "', '".join(token_matches)))
                else:
                    errors.append("'{}', a {}.{}, contains no identifiers matching the token pattern "
                        "of '{}' objects".format(expression,
                        model_class.__name__, attribute, ', '.join(object_type_names)))
                # break the analysis, because it is not clear how many tokens should be skipped
                break

            if 1 < len(id_matches):
                # like any typical lexer, pick the longest match
                id_matches_by_length = sorted(id_matches, key=lambda id_match: len(id_match.match_string))
                longest_length = len(id_matches_by_length[-1].match_string)
                longest_matches = set()
                while id_matches_by_length and len(id_matches_by_length[-1].match_string) == longest_length:
                    longest_matches.add(id_matches_by_length.pop())
                id_matches = longest_matches

            if 1 < len(id_matches):
                # error: multiple, maximal length model id_matches
                matches_error = ["'{}' as a {} id".format(id_val, model_type.__name__)
                    for model_type, _, id_val in sorted(id_matches, key=lambda id_match: id_match.model_type.__name__)]
                matches_error = ', '.join(matches_error)
                errors.append("'{}', a {}.{}, contains multiple model object id_matches: {}".format(expression,
                    model_class.__name__, attribute, matches_error))
                # again, break because it is not clear how many tokens should be skipped
                break

            else:
                # record an object ID token
                id_match = id_matches.pop()
                related_objects[id_match.model_type].append(id_match.match_string)
                wc_tokens.append(WcLangToken(TokCodes.wc_lang_obj_id, id_match.match_string, id_match.model_type))
                idx += len(id_match.token_pattern)
                continue

        # todo: perhaps ensure that the parsed expression can be eval'ed, assuming values in the range for related models
        # could even specify range in model declaration
        if errors:
            return (None, errors)
        else:
            # deterministically de-dupe all lists in related_objects
            for model_type, related_models in related_objects.items():
                related_objects[model_type] = det_dedupe(related_models)
            return (wc_tokens, related_objects)

    # the union of all valid_functions attributes for related Models
    valid_functions = (ceil, floor, exp, pow, log, log10, min, max)

    @staticmethod
    def eval_expr(dyn_model_obj, tokenized_expr, time, dynamic_model):
        """ Evaluate a Python expression in attribute `attribute` of object `obj`

        Called by the simulator when it calculates the value of a dynamic object, such as a
        DynamicObservable, or a RateLaw.

        Approach:
            * Replace references to related Models in `tokenized_expr` with their values
            * Join the elements of `tokenized_expr` into a Python expression
            * `eval` the Python expression

        Args:
            dyn_model_obj (:obj:`dyn_model_obj.Model`): a dynamic `wc_sim` `Model` instance whose expression is being evaluated
            tokenized_expr (:obj:`list` of `tuple`): the tokens in the deserialized expression
            time (:obj:`float`): the current simulation time
            dynamic_model (:obj:`wc_sim.DynamicModel`): a simulation's dynamical access method
            
        Returns:
            (:obj:`object`): the value of the expression at time `time`

        Raises:
            (:obj:`ValueError`): if the expression evaluation fails
        """
        # todo: ensure that all types of related Models can be evaluated through dynamic_model
        evaled_tokens = []
        for wc_token in tokenized_expr:
            if wc_token.tok_code != TokCodes.wc_lang_obj_id:
                evaled_tokens.append(wc_token.token_string)
            else:
                # evaluate the wc_lang_obj_id
                value = dynamic_model.eval_dynamic_obj(wc_token.model_type, wc_token.token_string, time)
                evaled_tokens.append(str(value))

        expression = ''.join(evaled_tokens)
        local_ns = {func.__name__: func for func in ExpressionUtils.valid_functions}
        # todo: what about k_cat and k_m?

        # get id whether it is a static attribute or a method, like in Species
        id = None
        if hasattr(dyn_model_obj, 'id'):
            id = getattr(dyn_model_obj, 'id')
            if callable(id):
                id = id()

        error_suffix = " cannot eval expression '{}' in {} with id {} at time {}; ".format(expression,
            dyn_model_obj.__class__.__name__, id, time)
        try:
            return eval(expression, {}, local_ns)
        except SyntaxError as error:
            raise ValueError("SyntaxError:" + error_suffix + str(error))
        except NameError as error:
            raise ValueError("NameError:" + error_suffix + str(error))
        except Exception as error:  # pragma: no cover
            raise ValueError("Exception:" + error_suffix + str(error))
