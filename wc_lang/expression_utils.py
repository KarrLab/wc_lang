""" Utilities for processing mathematical expressions used by wc_lang models

:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2017-10-23
:Copyright: 2016-2018, Karr Lab
:License: MIT
"""

import tokenize
import token
from six.moves import cStringIO
from io import BytesIO
from math import isnan, ceil, floor, exp, pow, log, log10
from collections import namedtuple

from wc_utils.util.list import det_dedupe
from wc_utils.util.enumerate import CaseInsensitiveEnum
from wc_utils.util.misc import DFSMAcceptor
import obj_model
import wc_lang

CONCENTRATIONS_DICT = 'concentrations'
PARAMETERS_DICT = 'parameters'

# TODO
'''
replace MakeModels
'''
'''
build
wc_lang:
    test WcLangExpression without accessing wc_lang.core
    right template for wc_sim_test
    use WcLangExpression to deserialize and validate all wc_lang expressions:
        last, convert ObjectiveFunction, and RateLawEquation (what about k_cat and k_m?)
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
    incorporate id into error_suffix in test_eval_expr
    test with real WC models
    perhaps make model_class a model, so that test_eval_expr can provide ids
    rename Function to Macro; need to change model files too -- argh
    expand Jupyter example
    support logical operators (or, and, not) in expressions, esp. StopConditions
cleanup:
    robust tests of run_results that examine the data
    remove transcode_and_check_rate_law_equations, which will be redundant (and ignores concentraion units)
    more docstrings and better naming for ExpressionVerifier
    use ExpressionMethods.make_expression_obj() wherever possible
    better error than "Not a linear expression of species and observables"
    do a better job of getting the plural in "if not used_model_type_attr.endswith('s'):"
    have valid_functions defined as sets, not tuples
    stop using & and remove RateLawUtils
    replace all validation and deserialization code
    fix test_find_shared_species
'''


class RateLawUtils(object):
    '''A set of static rate_law methods '''

    @staticmethod
    def transcode(rate_law_equation, species_ids, parameter_ids):
        '''Translate a `wc_lang.core.RateLawEquation` into a python expression that can be evaluated
        during a simulation.

        Args:
            rate_law_equation (:obj:`wc_lang.core.RateLawEquation`): a rate law equation
            species_ids (:obj:`set` of `str`): ids of the species that the rate law might use
            parameter_ids (:obj:`set` of `str`): ids of the parameters that the rate law might use

        Returns:
            The python expression, or None if the rate law doesn't have an equation

        Raises:
            :obj:`ValueError`: if `rate_law_equation` contains `__`, which increases its security risk, or
            if `rate_law_equation` refers to species not in `species` or parameters not in `parameters`
        '''
        def possible_specie_id(tokens):
            '''Indicate whether `tokens` begins with 4 tokens whose syntax id_matches a specie ID

            Specie IDs have the form `string[string]`, as documented in `wc_lang.core.Species`.

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
                :obj:`ValueError`: if `tokens` does not represent a specie in `species_ids`
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
                :obj:`ValueError`: if `tokens` does not represent a parameter in `parameter_ids`
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
    def eval_reaction_rate_laws(reaction, concentrations, parameters):
        '''Evaluate a reaction's rate laws at the given species concentrations and parameter values

        Args:
            reaction (:obj:`wc_lang.core.Reaction`): a Reaction instance
            concentrations (:obj:`dict` of :obj:`species_id` -> `float`):
                a dictionary of species concentrations
            parameters (:obj:`dict` of :obj:`parameter_id` -> `float`):
                a dictionary of parameter values

        Returns:
            (:obj:`list` of :obj:`float`): the reaction's forward and, perhaps, backward rates
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
            :obj:`ValueError`: if the rate law has a syntax error
            NameError: if the rate law references a specie whose concentration is not provided in
                `concentrations` or a parameter whose value is not provided in `parameters`
            Exception: if the rate law has other errors, such as a reference to an unknown function
        '''
        if not transcoded_equation:
            transcoded_equation = rate_law.equation.transcoded

        # todo: optimization: precompute local_ns for functions, k_cat, and k_m
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


class Error(Exception):
    """ Base class for exceptions in `WcLangExpression`

    Attributes:
        message (:obj:`str`): the exception's message
    """

    def __init__(self, message=None):
        super().__init__(message)


class WcLangExpressionError(Error):
    """ Exception raised for errors in `WcLangExpression`

    Attributes:
        message (:obj:`str`): the exception's message
    """

    def __init__(self, message=None):
        super().__init__(message)


class WcLangExpression(object):
    """ An expression in a wc_lang Model

    Expressions are currently (July, 2018) used in five `wc_lang` `Model`s: `RateLawEquation`, `Function`,
    `StopCondition` (which is just a special case of `Function` that returns a boolean), `ObjectiveFunction`,
    and `Observable`. These expressions are limited Python expressions with specific semantics:

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

    Attributes:
        model_class (:obj:`obj_model.Model`): the `wc_lang` `Model` which has an expression
        attribute (:obj:`str`): the attribute name of the expression in `model_class`
        expression (:obj:`str`): the expression defined in the wc_lang Model
        objects (:obj:`dict`): dict of wc_lang Models that might be referenced in expression; maps
            model type to a dict mapping ids to Model instances
        model_types (:obj:`dict`): wc_lang Model types that are needed by WcLangExpression code or can
            be used in expressions; a map from model name to wc_lang Model type
        tokens (:obj:`list` of :obj:`namedtuple`): a list of Python tokens generated by `tokenize.tokenize()`
        valid_used_models (:obj:`set`): wc_lang Models that `model_class` objects are allowed to use,
            as indicated in `model_class.Meta.valid_used_models`, intersected with `objects.keys()`
        valid_functions (:obj:`set`): the union of all `valid_functions` attributes for `objects`
        related_objects (:obj:`dict`): model instances that are referenced in `expression`; maps model type to
            dict that maps model id to model instance
        errors (:obj:`list` of :obj:`str`): errors found when parsing an `expression` fails
        wc_tokens (:obj:`list` of :obj:`WcLangToken`): tokens obtained when an `expression` is successfully
            `tokenize`d; if empty, then this `WcLangExpression` cannot use `eval_expr()`
    """

    # Function.identifier()
    FUN_TYPE_DISAMBIG_PATTTERN = (token.NAME, token.DOT, token.NAME, token.LPAR, token.RPAR)
    # ModelType.model_id
    MODEL_TYPE_DISAMBIG_PATTERN = (token.NAME, token.DOT, token.NAME)
    FUNCTION_PATTERN = (token.NAME, token.LPAR)

    # enumerate and detect Python tokens that are illegal in wc_lang expressions
    # TODO: consider the handful of other tokens that may also be illegal: COMMA, SEMI, TILDE, CIRCUMFLEX, and AT
    ILLEGAL_TOKEN_NAMES = ['ENDMARKER', 'NEWLINE', 'INDENT', 'DEDENT', 'COLON', 'LBRACE', 'RBRACE',
                           'PLUSEQUAL', 'MINEQUAL', 'STAREQUAL', 'SLASHEQUAL', 'PERCENTEQUAL', 'AMPEREQUAL',
                           'VBAREQUAL', 'CIRCUMFLEXEQUAL', 'LEFTSHIFTEQUAL', 'RIGHTSHIFTEQUAL', 'DOUBLESTAREQUAL',
                           'DOUBLESLASHEQUAL',
                           'ATEQUAL', 'RARROW', 'ELLIPSIS', 'AWAIT', 'ASYNC', 'ERRORTOKEN', 'N_TOKENS', 'NT_OFFSET']
    ILLEGAL_TOKENS = set()
    for illegal_name in ILLEGAL_TOKEN_NAMES:
        ILLEGAL_TOKENS.add(getattr(token, illegal_name))

    # These wc_lang models are needed by WcLangExpression, and must be provided in its constructor.
    REQUIRED_MODELS = ['Function']

    def __init__(self, model_class, attribute, expression, objects, given_model_types=None):
        """ Create an instance of WcLangExpression

        Args:
            model_class (:obj:`obj_model.Model`): the `wc_lang` `Model` which has an expression
            attribute (:obj:`str`): the attribute name of the expression in `model_class`
            expression (:obj:`str`): the expression defined in the wc_lang Model
            objects (:obj:`dict`): dict of wc_lang Models that might be referenced in expression; maps
                model type to a dict mapping ids to Model instances
            given_model_types (:obj:`list`): wc_lang Model types that are needed by WcLangExpression
                code or can be used in expressions

        Raises:
            :obj:`WcLangExpressionError`: if `model_class` is not a subclass of `obj_model.Model`,
                or lexical analysis of `expression` raises an exception,
                or `objects` includes model types that `model_class` should not reference
        """

        if not issubclass(model_class, obj_model.Model):
            raise WcLangExpressionError("model_class '{}' is not a subclass of obj_model.Model".format(
                model_class.__name__))
        if not hasattr(model_class.Meta, 'valid_used_models'):
            raise WcLangExpressionError("model_class '{}' doesn't have a 'Meta.valid_used_models' attribute".format(
                model_class.__name__))

        # Two categories of wc_lang model classes are needed, 1) the REQUIRED_MODELS used by
        # WcLangExpression code and 2) model types that can be referenced by model_class. Ensure that these
        # are provided in objects or given_model_types so their properties can be obtained.
        provided_model_types = set([model_type for model_type in objects.keys()])
        if given_model_types:
            provided_model_types = provided_model_types.union([model_type for model_type in given_model_types])
        needed_models = self.REQUIRED_MODELS + list(getattr(model_class.Meta, 'valid_used_models', []))
        missing_models = set(needed_models).difference([mt.__name__ for mt in provided_model_types])
        if missing_models:
            raise WcLangExpressionError(
                "model_class '{}': these needed models not in 'objects' or 'given_model_types': {}".format(
                model_class.__name__, missing_models))

        self.model_types = {model_type.__name__: model_type for model_type in provided_model_types}
        self.objects = objects

        self.valid_used_models = set()
        for valid_model_type_name in model_class.Meta.valid_used_models:
            valid_model_type = self.get_wc_lang_model_type(valid_model_type_name)
            self.valid_used_models.add(valid_model_type)
        for obj_type in self.valid_used_models:
            if not issubclass(obj_type, obj_model.Model):   # pragma    no cover
                raise WcLangExpressionError("objects entry '{}' is not a subclass of obj_model.Model".format(
                    obj_type.__name__))
        self.valid_functions = set()
        for obj_type in self.valid_used_models:
            if hasattr(obj_type.Meta, 'expression_model'):
                if hasattr(obj_type.Meta.expression_model.Meta, 'valid_functions'):
                    self.valid_functions.update(obj_type.Meta.expression_model.Meta.valid_functions)

        self.model_class = model_class
        self.attribute = attribute
        # strip leading and trailing whitespace from expression, which would create a bad token error
        self.expression = expression.strip()

        try:
            g = tokenize.tokenize(BytesIO(self.expression.encode('utf-8')).readline)
            # strip the leading ENCODING and trailing ENDMARKER tokens
            self.tokens = list(g)[1:-1]
        except tokenize.TokenError as e:
            raise WcLangExpressionError("parsing '{}', a {}.{}, creates a Python syntax error: '{}'".format(
                self.expression, self.model_class.__name__, self.attribute, str(e)))

        self.__reset_tokenization()

    def __reset_tokenization(self):
        """ Reset tokenization
        """
        self.related_objects = {}
        for model_type in self.valid_used_models:
            self.related_objects[model_type] = {}

        self.errors = []
        self.wc_tokens = []

    def get_wc_lang_model_type(self, model_type_name):
        """ Get the `wc_lang` model type corresponding to `model_type_name`

        Args:
            model_type_name (:obj:`str`): the name of a `wc_lang` model type

        Returns:
            :obj:`obj_model.Model`: the `wc_lang` model type with name `model_type_name`

        Raises:
            :obj:`WcLangExpressionError`: if no `wc_lang` model is named `model_type_name`
        """
        if model_type_name in self.model_types:
            return self.model_types[model_type_name]
        raise WcLangExpressionError("wc_lang model type '{}' not found".format(model_type_name))

    def get_used_model_type(self, model_type_name):
        """ Find the used `wc_lang` model type corresponding to `model_type_name`

        Args:
            model_type_name (:obj:`str`): the name of a purported `wc_lang` model type in an expression

        Returns:
            :obj:`object`: `None` if no model named `model_type_name` exists in `self.valid_used_models`,
                else the type of the model with that name
        """
        for model_type in self.valid_used_models:
            if model_type_name == model_type.__name__:
                return model_type
        return None

    def match_tokens(self, token_pattern, idx):
        """ Indicate whether `self.tokens` begins with a pattern of tokens that match `token_pattern`

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
                contains the casefold'ed value; identifier keys in `self.objects` must already be casefold'ed;
                default=False

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a disambugated id is found, return a `LexMatch` describing it.
        """
        fun_match = self.match_tokens(self.FUN_TYPE_DISAMBIG_PATTTERN, idx)
        if fun_match:
            possible_macro_id = self.tokens[idx+2].string
            if case_fold_match:
                possible_macro_id = possible_macro_id.casefold()
            # the disambiguation model type must be Function
            if self.tokens[idx].string != 'Function':
                return ("'{}', a {}.{}, contains '{}', which doesn't use 'Function' as a disambiguation "
                        "model type".format(self.expression, self.model_class.__name__, self.attribute, fun_match))
            # the identifier must be in the Function objects
            Function = self.get_wc_lang_model_type('Function')
            if Function not in self.valid_used_models or possible_macro_id not in self.objects[Function]:
                return "'{}', a {}.{}, contains '{}', which doesn't refer to a Function in 'objects'".format(
                    self.expression, self.model_class.__name__, self.attribute, fun_match)
            return LexMatch([WcLangToken(TokCodes.wc_lang_obj_id, fun_match, Function,
                                         possible_macro_id, self.objects[Function][possible_macro_id])],
                            len(self.FUN_TYPE_DISAMBIG_PATTTERN))

        disambig_model_match = self.match_tokens(self.MODEL_TYPE_DISAMBIG_PATTERN, idx)
        if disambig_model_match:
            disambig_model_type = self.tokens[idx].string
            possible_model_id = self.tokens[idx+2].string
            if case_fold_match:
                possible_model_id = possible_model_id.casefold()
            # the disambiguation model type cannot be Function
            if disambig_model_type == 'Function':
                return ("'{}', a {}.{}, contains '{}', which uses 'Function' as a disambiguation "
                        "model type but doesn't use Function syntax".format(self.expression, self.model_class.__name__,
                                                                            self.attribute, disambig_model_match))

            # the disambiguation model type must be in self.valid_used_models
            wc_lang_model_type = self.get_used_model_type(disambig_model_type)
            if wc_lang_model_type is None:
                return ("'{}', a {}.{}, contains '{}', but the disambiguation model type '{}' "
                        "cannot be referenced by '{}' expressions".format(
                            self.expression, self.model_class.__name__,
                            self.attribute, disambig_model_match, disambig_model_type,
                            self.model_class.__name__))

            if possible_model_id not in self.objects[wc_lang_model_type]:
                return "'{}', a {}.{}, contains '{}', but '{}' is not the id of a '{}'".format(
                    self.expression, self.model_class.__name__, self.attribute, disambig_model_match,
                    possible_model_id, disambig_model_type)

            return LexMatch([WcLangToken(TokCodes.wc_lang_obj_id, disambig_model_match, wc_lang_model_type,
                                         possible_model_id, self.objects[wc_lang_model_type][possible_model_id])],
                            len(self.MODEL_TYPE_DISAMBIG_PATTERN))

        # no match
        return None

    def related_object_id(self, idx, case_fold_match=False):
        """ Try to parse a related object `wc_lang` id from `self.tokens` at `idx`

        Different `wc_lang` objects match different Python token patterns. The default pattern
        is (token.NAME, ), but an object of type `model_type` can define a custom pattern in
        `model_type.Meta.token_pattern`, as Species does. Some patterns, like Species, may consume
        multiple Python tokens.

        Args:
            idx (:obj:`int`): current index into `tokens`
            case_fold_match (:obj:`bool`, optional): if set, casefold identifiers before matching;
                identifier keys in `self.objects` must already be casefold'ed; default=False

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a related object id is found, return a `LexMatch` describing it.
        """
        token_matches = set()
        id_matches = set()
        # find all matching related object ids
        for model_type in self.valid_used_models:
            token_pattern = (token.NAME, )
            if hasattr(model_type.Meta, 'token_pattern'):
                token_pattern = model_type.Meta.token_pattern
            match_string = self.match_tokens(token_pattern, idx)
            if match_string:
                token_matches.add(match_string)
                # skip model types that aren't in self.objects
                if model_type in self.objects:
                    # is match_string the ID of an instance in self.objects?
                    if case_fold_match:
                        if match_string.casefold() in self.objects[model_type]:
                            id_matches.add(IdMatch(model_type, token_pattern, match_string))
                    else:
                        if match_string in self.objects[model_type]:
                            id_matches.add(IdMatch(model_type, token_pattern, match_string))

        if not id_matches:
            if token_matches:
                return ("'{}', a {}.{}, contains the identifier(s) '{}', which aren't "
                        "the id(s) of an object in 'objects'".format(
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
                             self.objects[match.model_type][right_case_match_string])],
                len(match.token_pattern))

    def fun_call_id(self, idx, case_fold_match='unused'):
        """ Try to parse a Python math function call from `self.tokens` at `idx`

        Each `wc_lang` object `model_class` that contains an expression which can use Python math
        functions must define the set of allowed functions in `Meta.valid_functions` of the
        model_class Expression Model.

        Args:
            idx (:obj:`int`): current index into `self.tokens`
            case_fold_match (:obj:`str`, optional): ignored keyword; makes `WcLangExpression.tokenize()` simpler

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a function call is found, return a `LexMatch` describing it.
        """
        fun_match = self.match_tokens(self.FUNCTION_PATTERN, idx)
        if fun_match:
            fun_name = self.tokens[idx].string
            # FUNCTION_PATTERN is "identifier ("
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
                len(self.FUNCTION_PATTERN))

        # no match
        return None

    def tokenize(self, case_fold_match=False):
        """ Tokenize a Python expression in `self.expression`

        Args:
            case_fold_match (:obj:`bool`, optional): if set, casefold identifiers before matching;
                identifier keys in `self.objects` must already be casefold'ed; default = False

        Returns:
            (:obj:`tuple`): either `(None, :obj:list of :obj:str)` containing a list of errors, or
                `(:obj:list, :obj:dict)` containing a list of :obj:`WcLangToken`s and a dict of Model
                instances used by this list, grouped by Model type

        Raises:
            :obj:`WcLangExpressionError`: if `model_class` does not have a `Meta` attribute
        """
        self.__reset_tokenization()

        # detect and report bad tokens
        bad_tokens = set()
        for tok in self.tokens:
            if tok.exact_type in self.ILLEGAL_TOKENS:
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

        if self.errors:
            return (None, self.errors)
        return (self.wc_tokens, self.related_objects)

    def test_eval_expr(self, test_val=1.0):
        """ Test evaluate the Python expression in this `WcLangExpression`

        Called to validate this `WcLangExpression`.
        This expression must have been successfully `tokenize`d. All Models it uses are
        assumed to have the value `test_val`.

        Approach:
            * Replace references to used Models in `self.wc_tokens` with `test_val`
            * Join the elements of `self.wc_tokens` into a Python expression
            * `eval` the Python expression

        Args:
            test_val (:obj:`float`, optional): the value assumed for used Models

        Returns:
            (:obj:`object`): the value of the expression

        Raises:
            :obj:`WcLangExpressionError`: if the expression evaluation fails
        """
        # do not eval an expression that could not be tokenized
        if not self.wc_tokens:
            raise WcLangExpressionError("cannot evaluate '{}', as it not been successfully tokenized".format(
                self.expression))

        evaled_tokens = []
        idx = 0
        Function = self.get_wc_lang_model_type('Function')
        while idx < len(self.wc_tokens):
            wc_token = self.wc_tokens[idx]
            if wc_token.tok_code == TokCodes.wc_lang_obj_id:
                evaled_tokens.append(str(test_val))
                if wc_token.model_type == Function:
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
            raise WcLangExpressionError("SyntaxError:" + error_suffix + str(error))
        except NameError as error:  # pragma: no cover
            raise WcLangExpressionError("NameError:" + error_suffix + str(error))
        except Exception as error:  # pragma: no cover
            raise WcLangExpressionError("Exception:" + error_suffix + str(error))

    def __str__(self):
        rv = []
        rv.append("model_class: {}".format(self.model_class.__name__))
        rv.append("expression: '{}'".format(self.expression))
        rv.append("attribute: {}".format(self.attribute))
        rv.append("tokens: {}".format("'"+"', '".join([t.string for t in self.tokens])+"'"))
        rv.append("objects: {}".format(self.objects))
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
            return (False, "Not a linear expression of species and observables")


class LinearExpressionVerifier(ExpressionVerifier):
    """ Verify whether a sequence of tokens (`WcLangToken`s) describes a linear function of identifiers

    In particular, a valid linear expression must have the structure:
        * `(identifier | number '*' identifier) (('+' | '-') (identifier | number '*' identifier))*`
    """

    # Transitions in valid linear expression
    linear_expr_transitions = [   # (current state, message, next state)
        ('need number or id', (TokCodes.number, None), 'need *'),
        ('need *', (TokCodes.op, '*'), 'need id'),
        ('need id', (TokCodes.wc_lang_obj_id, None), 'need + | - | end'),
        ('need number or id', (TokCodes.wc_lang_obj_id, None), 'need + | - | end'),
        ('need + | - | end', (TokCodes.op, '+'), 'need number or id'),
        ('need + | - | end', (TokCodes.op, '-'), 'need number or id'),
        ('need + | - | end', (None, None), 'end')
    ]

    def __init__(self):
        super().__init__(start_state='need number or id', accepting_state='end',
                         transitions=self.linear_expr_transitions, empty_is_valid=True)

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
                return (False, "messages do not use token codes math_fun_id or other")
            if wc_lang_tok.tok_code == TokCodes.number:
                try:
                    float(wc_lang_tok.token_string)
                except ValueError as e:
                    return (False, str(e))
        return (True, True)

    def make_dfsa_messages(self, wc_lang_tokens):
        """ Convert a sequence of `WcLangToken`s into a list of messages for transitions in `linear_expr_transitions`

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
