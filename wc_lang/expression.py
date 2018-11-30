""" Utilities for processing mathematical expressions used by wc_lang models

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-30
:Copyright: 2016-2018, Karr Lab
:License: MIT
"""
import token
import tokenize
import obj_model
import wc_lang.core
from collections import namedtuple
from enum import Enum
from io import BytesIO
from wc_lang.core import InvalidObject, InvalidAttribute
from wc_lang.util import get_models
from wc_utils.util.misc import DFSMAcceptor


class WcTokenCodes(int, Enum):
    """ WcToken codes used in parsed expressions """
    wc_obj_id = 1
    math_func_id = 2
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
WcToken = namedtuple('WcToken', 'code, token_string, model_type, model_id, model')
# make model_type, model_id, and model optional: see https://stackoverflow.com/a/18348004
WcToken.__new__.__defaults__ = (None, None, None)
WcToken.__doc__ += ': WcToken in a parsed wc_lang expression'
WcToken.code.__doc__ = 'WcTokenCodes encoding'
WcToken.token_string.__doc__ = "The token's string"
WcToken.model_type.__doc__ = "When code is wc_obj_id, the wc_lang obj's type"
WcToken.model_id.__doc__ = "When code is wc_obj_id, the wc_lang obj's id"
WcToken.model.__doc__ = "When code is wc_obj_id, the wc_lang obj"


# result returned by a tokens lexer, like get_disambiguated_id()
LexMatch = namedtuple('LexMatch', 'wc_tokens, num_py_tokens')
LexMatch.__doc__ += ': result returned by a lexer method that matches a wc_lang expression element'
LexMatch.wc_tokens.__doc__ = 'List of WcLangTokens created'
LexMatch.num_py_tokens.__doc__ = 'Number of Python tokens consumed'


class Expression(object):
    """ Generic methods for mathematical expressions
    """

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return self.expression

    @classmethod
    def deserialize(cls, model_cls, value, objects):
        """ Deserialize expression

        Args:
            model_cls (:obj:`type`): expression class
            value (:obj:`str`): string representation of the mathematical expression, in a
                Python expression
            objects (:obj:`dict`): dictionary of objects which can be used in `expression`, grouped by model

        Returns:
            :obj:`tuple`: on error return (`None`, `InvalidAttribute`),
                otherwise return (object in this class with instantiated `_parsed_expression`, `None`)
        """
        # objects must contain all objects types in valid_models
        value = value or ''

        all_models = {model.__name__: model for model in get_models()}

        used_model_types = []
        for used_model in model_cls.Meta.valid_models:
            used_model_type = all_models[used_model]
            used_model_types.append(used_model_type)
        expr_field = 'expression'
        try:
            _parsed_expression = ParsedExpression(model_cls, expr_field, value, objects)
        except ParsedExpressionError as e:
            attr = model_cls.Meta.attributes['expression']
            return (None, InvalidAttribute(attr, [str(e)]))
        rv = _parsed_expression.tokenize()
        if rv[0] is None:
            attr = model_cls.Meta.attributes['expression']
            errors = rv[1]
            return (None, InvalidAttribute(attr, errors))
        _, used_objects = rv
        if model_cls not in objects:
            objects[model_cls] = {}
        if value in objects[model_cls]:
            obj = objects[model_cls][value]
        else:
            obj = model_cls(expression=value)
            objects[model_cls][value] = obj

            for attr_name, attr in model_cls.Meta.attributes.items():
                if isinstance(attr, obj_model.RelatedAttribute) and \
                        attr.related_class.__name__ in model_cls.Meta.valid_models:
                    attr_value = list(used_objects.get(attr.related_class, {}).values())
                    setattr(obj, attr_name, attr_value)
        obj._parsed_expression = _parsed_expression

        # check expression is linear
        obj._parsed_expression.is_linear, _ = LinearParsedExpressionVerifier().validate(
            obj._parsed_expression.wc_tokens)
        cls.set_lin_coeffs(obj)

        return (obj, None)

    @classmethod
    def set_lin_coeffs(cls, obj):
        """ Set the linear coefficients for the related objects

        Args:
            obj (:obj:`obj_model.Model`): expression object
        """
        model_cls = obj.__class__
        parsed_expr = obj._parsed_expression
        tokens = parsed_expr.wc_tokens
        is_linear = parsed_expr.is_linear

        if is_linear:
            default_val = 0.
        else:
            default_val = float('nan')

        parsed_expr.lin_coeffs = lin_coeffs = {}
        for attr_name, attr in model_cls.Meta.attributes.items():
            if isinstance(attr, obj_model.RelatedAttribute) and \
                    attr.related_class.__name__ in model_cls.Meta.valid_models:
                lin_coeffs[attr.related_class] = {}

        for related_class, related_objs in parsed_expr.related_objects.items():
            for related_obj in related_objs.values():
                lin_coeffs[related_class][related_obj] = default_val

        if not is_linear:
            return

        sense = 1.
        cur_coeff = 1.
        for token in tokens:
            if token.code == WcTokenCodes.op and token.token_string == '+':
                sense = 1.
                cur_coeff = 1.
            elif token.code == WcTokenCodes.op and token.token_string == '-':
                sense = -1.
                cur_coeff = 1.
            elif token.code == WcTokenCodes.number:
                cur_coeff = float(token.token_string)
            elif token.code == WcTokenCodes.wc_obj_id:
                lin_coeffs[token.model_type][token.model] += sense * cur_coeff

    @classmethod
    def validate(cls, model_obj, return_type=None, check_linear=False):
        """ Determine whether an expression model is valid by eval'ing its deserialized expression

        Args:
            model_obj (:obj:`Expression`): expression object
            return_type (:obj:`type`, optional): if provided, an expression's required return type
            check_linear (:obj:`bool`, optional): if :obj:`True`, validate that the expression is a
                linear function

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors in an `InvalidObject` instance
        """
        model_cls = model_obj.__class__

        # generate _parsed_expression
        objs = {}
        for related_attr_name, related_attr in model_cls.Meta.attributes.items():
            if isinstance(related_attr, obj_model.RelatedAttribute):
                objs[related_attr.related_class] = {
                    m.get_primary_attribute(): m for m in getattr(model_obj, related_attr_name)
                }
        try:
            model_obj._parsed_expression = ParsedExpression(model_obj.__class__, 'expression',
                                                            model_obj.expression, objs)
        except ParsedExpressionError as e:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, [str(e)])
            return InvalidObject(model_obj, [attr_err])

        is_valid, errors = model_obj._parsed_expression.tokenize()
        if is_valid is None:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, errors)
            return InvalidObject(model_obj, [attr_err])
        model_obj._parsed_expression.is_linear, _ = LinearParsedExpressionVerifier().validate(
            model_obj._parsed_expression.wc_tokens)
        cls.set_lin_coeffs(model_obj)

        # check related objects matches the tokens of the _parsed_expression
        related_objs = {}
        for related_attr_name, related_attr in model_cls.Meta.attributes.items():
            if isinstance(related_attr, obj_model.RelatedAttribute):
                related_model_objs = getattr(model_obj, related_attr_name)
                if related_model_objs:
                    related_objs[related_attr.related_class] = set(related_model_objs)

        token_objs = {}
        token_obj_ids = {}
        for token in model_obj._parsed_expression.wc_tokens:
            if token.model_type is not None:
                if token.model_type not in token_objs:
                    token_objs[token.model_type] = set()
                    token_obj_ids[token.model_type] = set()
                token_objs[token.model_type].add(token.model)
                token_obj_ids[token.model_type].add(token.token_string)

        if related_objs != token_objs:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, ['Related objects must match the tokens of the analyzed expression'])
            return InvalidObject(model_obj, [attr_err])

        # check expression is valid
        try:
            rv = model_obj._parsed_expression.test_eval()
            if return_type is not None:
                if not isinstance(rv, return_type):
                    attr = model_cls.Meta.attributes['expression']
                    attr_err = InvalidAttribute(attr,
                                                ["Evaluating '{}', a {} expression, should return a {} but it returns a {}".format(
                                                    model_obj.expression, model_obj.__class__.__name__,
                                                    return_type.__name__, type(rv).__name__)])
                    return InvalidObject(model_obj, [attr_err])
        except ParsedExpressionError as e:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, [str(e)])
            return InvalidObject(model_obj, [attr_err])

        # check expression is linear
        if check_linear and not model_obj._parsed_expression.is_linear:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, ['Expression must be linear'])
            return InvalidObject(model_obj, [attr_err])

        # return `None` to indicate valid object
        return None

    @staticmethod
    def make_expression_obj(model_type, expression, objs):
        """ Make an expression object

        Args:
            model_type (:obj:`type`): an `obj_model.Model` that uses a mathemetical expression, like
                `Function` and `Observable`
            expression (:obj:`str`): the expression used by the `model_type` being created
            objs (:obj:`dict` of `dict`): all objects that are referenced in `expression`

        Returns:
            :obj:`tuple`: if successful, (`obj_model.Model`, `None`) containing a new instance of
                `model_type`'s expression helper class; otherwise, (`None`, `InvalidAttribute`)
                reporting the error
        """
        expr_model_type = model_type.Meta.expression_model
        return expr_model_type.deserialize(expression, objs)

    @classmethod
    def make_obj(cls, model, model_type, id, expression, objs, allow_invalid_objects=False):
        """ Make a model that contains an expression by using its expression helper class

        For example, this uses `FunctionExpression` to make a `Function`.

        Args:
            model (:obj:`obj_model.Model`): a `wc_lang.core.Model` which is the root model
            model_type (:obj:`type`): an `obj_model.Model` that uses a mathemetical expression, like
                `Function` and `Observable`
            id (:obj:`str`): the id of the `model_type` being created
            expression (:obj:`str`): the expression used by the `model_type` being created
            objs (:obj:`dict` of `dict`): all objects that are referenced in `expression`
            allow_invalid_objects (:obj:`bool`, optional): if set, return object - not error - if
                the expression object does not validate

        Returns:
            :obj:`obj_model.Model` or `InvalidAttribute`: a new instance of `model_type`, or,
                if an error occurs, an `InvalidAttribute` reporting the error
        """
        expr_model_obj, error = cls.make_expression_obj(model_type, expression, objs)
        if error:
            return error
        error_or_none = expr_model_obj.validate()
        if error_or_none is not None and not allow_invalid_objects:
            return error_or_none
        related_name = model_type.Meta.attributes['model'].related_name
        related_in_model = getattr(model, related_name)
        new_obj = related_in_model.create(id=id, expression=expr_model_obj)
        return new_obj


class ParsedExpressionError(Exception):
    """ Exception raised for errors in `ParsedExpression`

    Attributes:
        message (:obj:`str`): the exception's message
    """

    def __init__(self, message=None):
        """
        Args:
            message (:obj:`str`, optional): the exception's message
        """
        super().__init__(message)


# enumerate and detect Python tokens that are illegal in wc_lang expressions
ILLEGAL_TOKENS_NAMES = ('ENDMARKER', 'NEWLINE', 'INDENT', 'DEDENT', 'COLON', 'LBRACE', 'RBRACE',
                        'PLUSEQUAL', 'MINEQUAL', 'STAREQUAL', 'SLASHEQUAL', 'PERCENTEQUAL', 'AMPEREQUAL', 'VBAREQUAL',
                        'CIRCUMFLEXEQUAL', 'LEFTSHIFTEQUAL', 'RIGHTSHIFTEQUAL', 'DOUBLESTAREQUAL', 'DOUBLESLASHEQUAL',
                        'ATEQUAL', 'RARROW', 'ELLIPSIS', 'AWAIT', 'ASYNC', 'ERRORTOKEN', 'N_TOKENS', 'NT_OFFSET',
                        'SEMI', 'TILDE', 'AT')
ILLEGAL_TOKENS = set()
for illegal_name in ILLEGAL_TOKENS_NAMES:
    ILLEGAL_TOKENS.add(getattr(token, illegal_name))


class ParsedExpression(object):
    """ An expression in a wc_lang Model

    Expressions are currently (July, 2018) used in five `wc_lang` `Model`s: `RateLawExpression`, `Function`,
    `StopCondition` (which is just a special case of `Function` that returns a boolean), `DfbaObjective`,
    and `Observable`. These expressions are limited Python expressions with specific semantics:

    * They must be syntactically correct Python.
    * No Python keywords, strings, or tokens that do not belong in expressions are allowed.
    * All Python identifiers must be the ID of an object in a whole-cell model, or components of
        an ID of an object in the model, or the name of a function in the `math` package. Objects in the model
        are provided in `_objs`, and the allowed subset of functions in `math` must be provided in an
        iterator in the `valid_functions` attribute of the `Meta` class of a model whose whose expression
        is being processed.
    * Currently (July, 2018), identifiers may refer to `Species`s, `Parameter`s, `Observable`s, `Reaction`s,
        `Observable`'s and `DfbaNetReaction`s.
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
        * `DfbaNetReaction`: its current flux
    * The modeller is responsible for ensuring that units in expressions are internally consistent and appropriate
        for the expression's use

    Attributes:
        model_cls (:obj:`type`): the `wc_lang` `Model` which has an expression
        attr (:obj:`str`): the attribute name of the expression in `model_cls`
        expression (:obj:`str`): the expression defined in the wc_lang Model
        py_tokens (:obj:`list` of :obj:`namedtuple`): a list of Python tokens generated by `tokenize.tokenize()`
        _objs (:obj:`dict`): dict of wc_lang Models that might be referenced in expression; maps
            model type to a dict mapping ids to Model instances
        valid_models (:obj:`set`): wc_lang Models that `model_cls` objects are allowed to use,
            as indicated in `model_cls.Meta.valid_models`, intersected with `_objs.keys()`
            might be referenced in expression; maps
        valid_functions (:obj:`set`): the union of all `valid_functions` attributes for `_objs`
        related_objects (:obj:`dict`): models that are referenced in `expression`; maps model type to
            dict that maps model id to model instance
        lin_coeffs (:obj:`dict`): linear coefficients of models that are referenced in `expression`;
            maps model type to dict that maps models to coefficients
        errors (:obj:`list` of :obj:`str`): errors found when parsing an `expression` fails
        wc_tokens (:obj:`list` of :obj:`WcToken`): tokens obtained when an `expression` is successfully
            `tokenize`d; if empty, then this `ParsedExpression` cannot use `eval()`
    """

    # Function.identifier()
    FUNC_TYPE_DISAMBIG_PATTERN = (token.NAME, token.DOT, token.NAME, token.LPAR, token.RPAR)
    # ModelType.model_id
    MODEL_TYPE_DISAMBIG_PATTERN = (token.NAME, token.DOT, token.NAME)
    FUNC_PATTERN = (token.NAME, token.LPAR)

    def __init__(self, model_cls, attr, expression, objs):
        """ Create an instance of ParsedExpression

        Args:
            model_cls (:obj:`type`): the `wc_lang` `Model` which has an expression
            attr (:obj:`str`): the attribute name of the expression in `model_cls`
            expression (:obj:`str`): the expression defined in the wc_lang Model
            objs (:obj:`dict`): dictionary of model objects (instances of :obj:`obj_model.Model`) organized
                by their type

        Raises:
            :obj:`ParsedExpressionError`: if `model_cls` is not a subclass of `obj_model.Model`,
                or lexical analysis of `expression` raises an exception,
                or `objs` includes model types that `model_cls` should not reference
        """
        if not issubclass(model_cls, obj_model.Model):
            raise ParsedExpressionError("model_cls '{}' is not a subclass of obj_model.Model".format(
                model_cls.__name__))
        if not hasattr(model_cls.Meta, 'valid_models'):
            raise ParsedExpressionError("model_cls '{}' doesn't have a 'Meta.valid_models' attribute".format(
                model_cls.__name__))
        self.valid_models = set()
        for valid_model_type_name in model_cls.Meta.valid_models:
            valid_model_type = getattr(wc_lang.core, valid_model_type_name)
            if valid_model_type in objs:
                self.valid_models.add(valid_model_type)
        for obj_type in self.valid_models:
            if not issubclass(obj_type, obj_model.Model):   # pragma    no cover
                raise ParsedExpressionError("objs entry '{}' is not a subclass of obj_model.Model".format(
                    obj_type.__name__))
        self.valid_functions = set()
        if hasattr(model_cls.Meta, 'valid_functions'):
            self.valid_functions.update(model_cls.Meta.valid_functions)

        self._objs = objs
        self.model_cls = model_cls
        self.attr = attr
        # strip leading and trailing whitespace from expression, which would create a bad token error
        self.expression = expression.strip()

        try:
            g = tokenize.tokenize(BytesIO(self.expression.encode('utf-8')).readline)
            # strip the leading ENCODING and trailing ENDMARKER tokens
            self.py_tokens = list(g)[1:-1]
        except tokenize.TokenError as e:
            raise ParsedExpressionError("parsing '{}', a {}.{}, creates a Python syntax error: '{}'".format(
                self.expression, self.model_cls.__name__, self.attr, str(e)))

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

    def get_model_type(self, name):
        """ Find the `wc_lang` model type corresponding to `name`

        Args:
            name (:obj:`str`): the name of a purported `wc_lang` model type in an expression

        Returns:
            :obj:`object`: `None` if no model named `name` exists in `self.valid_models`,
                else the type of the model with that name
        """
        for model_type in self.valid_models:
            if name == model_type.__name__:
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
        if len(self.py_tokens)-idx < len(token_pattern):
            return False
        for tok_idx, token_pat_num in enumerate(token_pattern):
            if self.py_tokens[idx+tok_idx].exact_type != token_pat_num:
                return False
            # because a wc_lang ID shouldn't contain white space, do not allow it between the self.py_tokens
            # that match token_pattern
            if 0 < tok_idx and self.py_tokens[idx+tok_idx-1].end != self.py_tokens[idx+tok_idx].start:
                return False
        match_val = ''.join([self.py_tokens[idx+i].string for i in range(len(token_pattern))])
        return match_val

    def get_disambiguated_id(self, idx, case_fold_match=False):
        """ Try to parse a disambiguated `wc_lang` id from `self.py_tokens` at `idx`

        Look for a disambugated id (either a Function written as `Function.identifier()`, or a
        Model written as `ModelType.model_id`). If tokens do not match, return `None`. If tokens match,
        but their values are wrong, return an error `str`.
        If a disambugated id is found, return a `LexMatch` describing it.

        Args:
            idx (:obj:`int`): current index into `tokens`
            case_fold_match (:obj:`bool`, optional): if set, `casefold()` identifiers before matching;
                in a `WcToken`, `token_string` retains the original expression text, while `model_id`
                contains the casefold'ed value; identifier keys in `self._objs` must already be casefold'ed;
                default=False

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a disambugated id is found, return a `LexMatch` describing it.
        """
        func_match = self.match_tokens(self.FUNC_TYPE_DISAMBIG_PATTERN, idx)
        if func_match:
            possible_macro_id = self.py_tokens[idx+2].string
            if case_fold_match:
                possible_macro_id = possible_macro_id.casefold()
            # the disambiguation model type must be Function
            if self.py_tokens[idx].string != wc_lang.core.Function.__name__:
                return ("'{}', a {}.{}, contains '{}', which doesn't use 'Function' as a disambiguation "
                        "model type".format(self.expression, self.model_cls.__name__, self.attr, func_match))
            # the identifier must be in the Function objects
            if wc_lang.core.Function not in self.valid_models or possible_macro_id not in self._objs[wc_lang.core.Function]:
                return "'{}', a {}.{}, contains '{}', which doesn't refer to a Function".format(
                    self.expression, self.model_cls.__name__, self.attr, func_match)
            return LexMatch([WcToken(WcTokenCodes.wc_obj_id, func_match, wc_lang.core.Function,
                                     possible_macro_id, self._objs[wc_lang.core.Function][possible_macro_id])],
                            len(self.FUNC_TYPE_DISAMBIG_PATTERN))

        disambig_model_match = self.match_tokens(self.MODEL_TYPE_DISAMBIG_PATTERN, idx)
        if disambig_model_match:
            disambig_model_type = self.py_tokens[idx].string
            possible_model_id = self.py_tokens[idx+2].string
            if case_fold_match:
                possible_model_id = possible_model_id.casefold()
            # the disambiguation model type cannot be Function
            if disambig_model_type == wc_lang.core.Function.__name__:
                return ("'{}', a {}.{}, contains '{}', which uses 'Function' as a disambiguation "
                        "model type but doesn't use Function syntax".format(self.expression, self.model_cls.__name__,
                                                                            self.attr, disambig_model_match))

            # the disambiguation model type must be in self.valid_models
            model_type = self.get_model_type(disambig_model_type)
            if model_type is None:
                return ("'{}', a {}.{}, contains '{}', but the disambiguation model type '{}' "
                        "cannot be referenced by '{}' expressions".format(
                            self.expression, self.model_cls.__name__,
                            self.attr, disambig_model_match, disambig_model_type,
                            self.model_cls.__name__))

            if possible_model_id not in self._objs[model_type]:
                return "'{}', a {}.{}, contains '{}', but '{}' is not the id of a '{}'".format(
                    self.expression, self.model_cls.__name__, self.attr, disambig_model_match,
                    possible_model_id, disambig_model_type)

            return LexMatch([WcToken(WcTokenCodes.wc_obj_id, disambig_model_match, model_type,
                                     possible_model_id, self._objs[model_type][possible_model_id])],
                            len(self.MODEL_TYPE_DISAMBIG_PATTERN))

        # no match
        return None

    def get_related_obj_id(self, idx, case_fold_match=False):
        """ Try to parse a related object `wc_lang` id from `self.py_tokens` at `idx`

        Different `wc_lang` objects match different Python token patterns. The default pattern
        is (token.NAME, ), but an object of type `model_type` can define a custom pattern in
        `model_type.Meta.token_pattern`, as Species does. Some patterns may consume multiple Python tokens.

        Args:
            idx (:obj:`int`): current index into `py_tokens`
            case_fold_match (:obj:`bool`, optional): if set, casefold identifiers before matching;
                identifier keys in `self._objs` must already be casefold'ed; default=False

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
                    if match_string.casefold() in self._objs[model_type]:
                        id_matches.add(IdMatch(model_type, token_pattern, match_string))
                else:
                    if match_string in self._objs[model_type]:
                        id_matches.add(IdMatch(model_type, token_pattern, match_string))

        if not id_matches:
            if token_matches:
                return ("'{}', a {}.{}, contains the identifier(s) '{}', which aren't "
                        "the id(s) of an object".format(
                            self.expression, self.model_cls.__name__,
                            self.attr, "', '".join(token_matches)))
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
                self.expression, self.model_cls.__name__, self.attr, matches_error)

        else:
            # return a lexical match about a related id
            match = id_matches.pop()
            right_case_match_string = match.match_string
            if case_fold_match:
                right_case_match_string = match.match_string.casefold()
            return LexMatch(
                [WcToken(WcTokenCodes.wc_obj_id, match.match_string, match.model_type, right_case_match_string,
                         self._objs[match.model_type][right_case_match_string])],
                len(match.token_pattern))

    def get_func_call_id(self, idx, case_fold_match='unused'):
        """ Try to parse a Python math function call from `self.py_tokens` at `idx`

        Each `wc_lang` object `model_cls` that contains an expression which can use Python math
        functions must define the set of allowed functions in `Meta.valid_functions` of the
        model_cls Expression Model.

        Args:
            idx (:obj:`int`): current index into `self.py_tokens`
            case_fold_match (:obj:`str`, optional): ignored keyword; makes `ParsedExpression.tokenize()` simpler

        Returns:
            :obj:`object`: If tokens do not match, return `None`. If tokens match,
                but their values are wrong, return an error `str`.
                If a function call is found, return a `LexMatch` describing it.
        """
        func_match = self.match_tokens(self.FUNC_PATTERN, idx)
        if func_match:
            func_name = self.py_tokens[idx].string
            # FUNC_PATTERN is "identifier ("
            # the closing paren ")" will simply be encoded as a WcToken with code == op

            # are Python math functions defined?
            if not hasattr(self.model_cls.Meta, 'valid_functions'):
                return ("'{}', a {}.{}, contains the func name '{}', but {}.Meta doesn't "
                        "define 'valid_functions'".format(self.expression,
                                                          self.model_cls.__name__, self.attr, func_name, self.model_cls.__name__))

            function_ids = set([f.__name__ for f in self.model_cls.Meta.valid_functions])

            # is the function allowed?
            if func_name not in function_ids:
                return ("'{}', a {}.{}, contains the func name '{}', but it isn't in "
                        "{}.Meta.valid_functions: {}".format(self.expression, self.model_cls.__name__,
                                                             self.attr, func_name, self.model_cls.__name__, ', '.join(function_ids)))

            # return a lexical match about a math function
            return LexMatch(
                [WcToken(WcTokenCodes.math_func_id, func_name), WcToken(WcTokenCodes.op, '(')],
                len(self.FUNC_PATTERN))

        # no match
        return None

    def tokenize(self, case_fold_match=False):
        """ Tokenize a Python expression in `self.expression`

        Args:
            case_fold_match (:obj:`bool`, optional): if set, casefold identifiers before matching;
                identifier keys in `self._objs` must already be casefold'ed; default = False

        Returns:
            (:obj:`tuple`): either `(:obj:`None`, :obj:`list` of :obj:`str)` containing a list of errors, or
                `(:obj:`list`, :obj:`dict`)` containing a list of :obj:`WcToken`s and a dict of Model
                instances used by this list, grouped by Model type

        Raises:
            :obj:`ParsedExpressionError`: if `model_cls` does not have a `Meta` attribute
        """
        self.__reset_tokenization()

        # detect and report bad tokens
        bad_tokens = set()
        for tok in self.py_tokens:
            print(tok.exact_type, tok.string)
            if tok.exact_type in ILLEGAL_TOKENS:
                if tok.string and tok.string != ' ':
                    bad_tokens.add(tok.string)
                else:
                    bad_tokens.add(token.tok_name[tok.type])
        if bad_tokens:
            self.errors.append("'{}', a {}.{}, contains bad token(s): '{}'".format(
                self.expression, self.model_cls.__name__,
                self.attr, "', '".join(bad_tokens)))
            return (None, self.errors)

        idx = 0
        while idx < len(self.py_tokens):

            # categorize token codes
            wc_token_code = WcTokenCodes.other
            if self.py_tokens[idx].type == token.OP:
                wc_token_code = WcTokenCodes.op
            elif self.py_tokens[idx].type == token.NUMBER:
                wc_token_code = WcTokenCodes.number

            # a token that isn't an identifier needs no processing
            if self.py_tokens[idx].type != token.NAME:
                # record non-identifier token
                self.wc_tokens.append(WcToken(wc_token_code, self.py_tokens[idx].string))
                idx += 1
                continue

            matches = []
            tmp_errors = []
            for get_wc_lex_el in [self.get_related_obj_id, self.get_disambiguated_id, self.get_func_call_id]:
                result = get_wc_lex_el(idx, case_fold_match=case_fold_match)
                if result is not None:
                    if isinstance(result, str):
                        tmp_errors.append(result)
                    elif isinstance(result, LexMatch):
                        matches.append(result)
                    else:   # pragma no cover
                        raise ParsedExpressionError("Result is neither str nor LexMatch '{}'".format(result))

            # should find either matches or errors
            if not (matches or tmp_errors):
                raise ParsedExpressionError("No matches or errors found in '{}'".format(self.expression))
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
            if len(longest_matches) > 1:
                raise ParsedExpressionError("Multiple longest matches: '{}'".format(longest_matches))

            # good match
            # advance idx to the next token
            # record match data in self.wc_tokens and self.related_objects
            match = longest_matches.pop()
            idx += match.num_py_tokens
            wc_tokens = match.wc_tokens
            self.wc_tokens.extend(wc_tokens)
            for wc_token in wc_tokens:
                if wc_token.code == WcTokenCodes.wc_obj_id:
                    self.related_objects[wc_token.model_type][wc_token.model_id] = wc_token.model

        # detect ambiguous tokens
        valid_function_names = [func.__name__ for func in self.valid_functions]
        for wc_token in self.wc_tokens:
            if wc_token.code in [WcTokenCodes.wc_obj_id, WcTokenCodes.math_func_id]:
                matching_items = []

                for model_type in self.valid_models:
                    if wc_token.token_string in self._objs[model_type]:
                        matching_items.append(model_type.__name__)

                if wc_token.token_string in valid_function_names:
                    matching_items.append('function')

                if len(matching_items) > 1:
                    self.errors.append('WcToken `{}` is ambiguous. WcToken matches a {} and a {}.'.format(
                        wc_token.token_string, ', a '.join(matching_items[0:-1]), matching_items[-1]))

        if self.errors:
            return (None, self.errors)
        return (self.wc_tokens, self.related_objects)

    def test_eval(self, val=1.0):
        """ Test evaluate the expression with the value of all used models equal to `test_val`.

        Called to validate this :obj:`ParsedExpression`.

        Args:
            val (:obj:`float`, optional): the value assumed for used models

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
                model_vals[model_type][model] = val

        return self.eval(model_vals)

    def eval(self, vals):
        """ Evaluate the expression

        This expression must have been successfully `tokenize`d.

        Approach:
            * Replace references to used Models in `self.wc_tokens` with values
            * Join the elements of `self.wc_tokens` into a Python expression
            * `eval` the Python expression

        Args:
            vals (:obj:`dict` of `Model): dictionary of values of the models; maps
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
            if wc_token.code == WcTokenCodes.wc_obj_id:
                evaled_tokens.append(str(vals[wc_token.model_type][wc_token.model]))
                if wc_token.model_type == wc_lang.core.Function:
                    # skip past the following ( ) tokens -- they're just syntactic sugar for Functions
                    idx += 2
            else:
                evaled_tokens.append(wc_token.token_string)
            idx += 1

        expression = ' '.join(evaled_tokens)
        local_ns = {func.__name__: func for func in self.valid_functions}

        error_suffix = " cannot eval expression '{}' in {}; ".format(expression,
                                                                     self.model_cls.__name__)

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
        rv.append("model_cls: {}".format(self.model_cls.__name__))
        rv.append("expression: '{}'".format(self.expression))
        rv.append("attr: {}".format(self.attr))
        rv.append("py_tokens: {}".format("'"+"', '".join([t.string for t in self.py_tokens])+"'"))
        rv.append("related_objects: {}".format(self.related_objects))
        rv.append("errors: {}".format(self.errors))
        rv.append("wc_tokens: {}".format(self.wc_tokens))
        return '\n'.join(rv)


class ParsedExpressionVerifier(object):
    """ Verify whether a sequence of `WcToken` tokens

    An `ParsedExpressionVerifier` consists of two parts:

    * An optional method `validate_tokens` that examines the content of individual tokens
      and returns `(True, True)` if they are all valid, or (`False`, error) otherwise. It can be
      overridden by subclasses.
    * A `DFSMAcceptor` which determines whether the tokens describe a particular pattern
      `validate()` combines these parts.

    Attributes:
        dfsm_acceptor (:obj:`DFSMAcceptor`): the DFSM acceptor
        empty_is_valid (:obj:`bool`): if set, then an empty sequence of tokens is valid
    """

    def __init__(self, start_state, accepting_state, transitions, empty_is_valid=False):
        """
        Args:
            start_state (:obj:`object`): a DFSM's start state
            accepting_state (:obj:`object`): a DFSM must be in this state to accept a message sequence
            transitions (:obj:`iterator` of `tuple`): transitions, an iterator of
                (state, message, next state) tuples
            empty_is_valid (:obj:`bool`, optional): if set, then an empty sequence of tokens is valid
        """
        self.dfsm_acceptor = DFSMAcceptor(start_state, accepting_state, transitions)
        self.empty_is_valid = empty_is_valid

    def validate_tokens(self, tokens):
        """ Check whether the content of a sequence of `WcToken`s is valid

        In particular, all numbers in `tokens` must be floats, and all token codes must not
        be `math_func_id` or `other`.

        Args:
            tokens (:obj:`iterator` of `WcToken`): sequence of `WcToken`s

        Returns:
            :obj:`tuple`: (`False`, error) if `tokens` cannot be a linear expression, or
                (`True`, `True`) if it can
        """
        return (True, True)

    def make_dfsa_messages(self, tokens):
        """ Convert a sequence of `WcToken`s into a list of messages for transitions

        Args:
            tokens (:obj:`iterator` of `WcToken`): sequence of `WcToken`s

        Returns:
            :obj:`object`: `None` if `tokens` cannot be converted into a sequence of messages,
                or a `list` of `tuple` of pairs (token code, message modifier)
        """
        messages = []
        for token in tokens:
            messages.append((token.code, None))
        return messages

    def validate(self, tokens):
        """ Indicate whether `tokens` is valid

        Args:
            tokens (:obj:`iterator` of `WcToken`): sequence of `WcToken`s

        Returns:
            :obj:`tuple`: (`False`, error) if `tokens` is valid, or (`True`, `None`) if it is
        """
        if self.empty_is_valid and not tokens:
            return (True, None)
        valid, error = self.validate_tokens(tokens)
        if not valid:
            return (False, error)
        dfsa_messages = self.make_dfsa_messages(tokens)
        if DFSMAcceptor.ACCEPT == self.dfsm_acceptor.run(dfsa_messages):
            return (True, None)
        else:
            return (False, "Not a linear expression")


class LinearParsedExpressionVerifier(ParsedExpressionVerifier):
    """ Verify whether a sequence of tokens (`WcToken`s) describes a linear function of identifiers

    In particular, a valid linear expression must have the structure:
        * `(identifier | number '*' identifier) (('+' | '-') (identifier | number '*' identifier))*`
    """

    # Transitions in valid linear expression
    TRANSITIONS = [   # (current state, message, next state)
        ('need number or id', (WcTokenCodes.number, None), 'need * id'),
        ('need * id', (WcTokenCodes.op, '*'), 'need id'),
        ('need id', (WcTokenCodes.wc_obj_id, None), 'need + | - | end'),
        ('need number or id', (WcTokenCodes.wc_obj_id, None), 'need + | - | end'),
        ('need + | - | end', (WcTokenCodes.op, '+'), 'need number or id'),
        ('need + | - | end', (WcTokenCodes.op, '-'), 'need number or id'),
        ('need + | - | end', (None, None), 'end'),
    ]

    def __init__(self):
        super().__init__(start_state='need number or id', accepting_state='end',
                         transitions=self.TRANSITIONS, empty_is_valid=True)

    def validate_tokens(self, tokens):
        """ Check whether the content of a sequence of `WcToken`s is valid

        In particular, all numbers in `tokens` must be floats, and all token codes must not
        be `math_func_id` or `other`.

        Args:
            tokens (:obj:`iterator` of `WcToken`): sequence of `WcToken`s

        Returns:
            :obj:`tuple`: (`False`, error) if `tokens` cannot be a linear expression, or
                (`True`, `True`) if it can
        """
        for token in tokens:
            if token.code in set([WcTokenCodes.math_func_id, WcTokenCodes.other]):
                return (False, "messages do not use token codes `math_func_id` or `other`")
            if token.code == WcTokenCodes.number:
                try:
                    float(token.token_string)
                except ValueError as e:
                    return (False, str(e))

        return (True, True)

    def make_dfsa_messages(self, tokens):
        """ Convert a sequence of `WcToken`s into a list of messages for transitions in
        :obj:`LinearParsedExpressionVerifier.TRANSITIONS`

        Args:
            tokens (:obj:`iterator` of `WcToken`): sequence of `WcToken`s

        Returns:
            :obj:`object`: :obj:`None` if `tokens` cannot be converted into a sequence of messages
                to validate a linear expression, or a :obj:`list` of :obj:`tuple` of pairs (token code, message modifier)
        """
        messages = []
        for token in tokens:
            message_tok_code = token.code
            if token.code == WcTokenCodes.wc_obj_id:
                message_modifier = None
            elif token.code == WcTokenCodes.number:
                message_modifier = None
            elif token.code == WcTokenCodes.op:
                message_modifier = token.token_string
            else:
                return None
            messages.append((message_tok_code, message_modifier))
        messages.append((None, None))
        return messages
