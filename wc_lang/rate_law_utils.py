'''Static methods which handle rate laws and their expressions

:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2017-10-23
:Copyright: 2016-2017, Karr Lab
:License: MIT
'''
import tokenize
import token
from six.moves import cStringIO
from math import isnan

import wc_lang

CONCENTRATIONS_DICT = 'concentrations'


# TODO(Arthur): use this to deserialize rate law expressions
# TODO(Arthur): add functionality to create rate law's modifiers attribute
class RateLawUtils(object):
    '''A set of static rate_law methods '''

    @staticmethod
    def transcode(rate_law_equation, species):
        '''Translate a `wc_lang.core.RateLawEquation` into a python expression that can be evaluated
        during a simulation.

        Args:
            rate_law_equation (:obj:`wc_lang.core.RateLawEquation`): a rate law equation
            species (:obj:`set` of `wc_lang.core.Species`): the species that use the rate law

        Returns:
            The python expression, or None if the rate law doesn't have an equation

        Raises:
            ValueError: If `rate_law_equation` contains `__`, which increases its security risk, or
            if `rate_law_equation` refers to species not in `species`
        '''
        def possible_specie_id(tokens):
            '''Indicate whether `tokens` begins with 4 tokens whose syntax matches a specie ID

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

        if '__' in rate_law_equation.expression:
            raise ValueError("Security risk: rate law expression '{}' contains '__'.".format(
                rate_law_equation.expression))

        rate_law_expression = rate_law_equation.expression
        species_ids = set([specie.id() for specie in species])

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
                                                                              submodel.get_species())

    @staticmethod
    def eval_reaction_rate_laws(reaction, concentrations):
        '''Evaluate a reaction's rate laws at the given species concentrations

        Args:
            reaction (:obj:`wc_lang.core.Reaction`): a Reaction instance
            concentrations (:obj:`dict` of :obj:`species_id` -> `float`):
                a dictionary of species concentrations

        Returns:
            (:obj:`list` of `float`): the reaction's forward and, perhaps, backward rates

        Raises:
            ValueError: if the reaction's rate law has a syntax error
            NameError: if the rate law references a specie whose concentration is not provided in
                `concentrations`
            Exception: if the rate law has other errors, such as a reference to an unknown function
        '''
        rates = []
        for rate_law in reaction.rate_laws:
            rates.append(RateLawUtils.eval_rate_law(rate_law, concentrations))
        return rates

    @staticmethod
    def eval_rate_law(rate_law, concentrations, transcoded_equation=None):
        '''Evaluate a rate law at the given species concentrations

        Args:
            rate_law (:obj:`wc_lang.core.RateLaw`): a RateLaw instance
            concentrations (:obj:`dict` of :obj:`species_id` -> `float`):
                a dictionary of species concentrations
            transcoded_equation (:obj:`str`, optional): the rate law's transcoded_equation; if not
                provided, will be taken from rate_law.equation.transcoded

        Returns:
            (:obj:`float`): the rate law's rate

        Raises:
            ValueError: if the rate law has a syntax error
            NameError: if the rate law references a specie whose concentration is not provided in
                `concentrations`
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
