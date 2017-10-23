'''Static methods which handle rate laws and their expressions

:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Date: 2017-10-23
:Copyright: 2016-2017, Karr Lab
:License: MIT
'''
import tokenize, token
from six.moves import cStringIO
from math import isnan

from wc_lang.core import RateLaw, RateLawEquation, Species, Reaction, Model

CONCENTRATIONS_DICT = 'concentrations'

class RateLawUtils(object):
    '''A set of static rate_law methods '''

    @staticmethod
    def transcode(rate_law, species):
        '''Translate a `RateLaw` into a python expression that can be evaluated
        during a simulation.

        Args:
            rate_law (:obj:`RateLaw`): a rate law
            species (:obj:`set` of `Species`): the species that use the rate law

        Returns:
            The python expression, or None if the rate law doesn't have an equation

        Raises:
            ValueError: If `rate_law` contains `__`, which increases its security risk, or
            if `rate_law` refers to species not in `species`
        '''
        def possible_specie_id(tokens):
            '''Determine whether `tokens` contains a specie id of the form 'string[string]'

            See `Species.serialize()` for the format of specie ids.

            Args:
                tokens (:obj:`list` of (token_num, token_val)): a list of Python tokens

            Returns:
                True if `tokens` might contain a specie id
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
            '''Translate `tokens` into a python expression that can be evaluated during a simulation

            Args:
                tokens (:obj:`list` of (token_num, token_val)): a list of 4 Python tokens that
                    should comprise a specie id
                species_ids (:obj:`set`): ids of the species used by the rate law expression
                rate_law_expression (:obj:`string`): the rate law expression being transcoded

            Returns:
                (:obj:`string`): a Python expression, transcoded to look up the specie concentration
                    in `concentrations[]`

            Raises:
                ValueError: if `tokens` does not represent a specie in `species_ids`
            '''
            tokvals = [token_tmp[1] for token_tmp in tokens]
            # TODO: centralize the format species ids in a pair of serialize & parse utilities
            parsed_id = "{}[{}]".format(tokvals[0], tokvals[2])
            if parsed_id in species_ids:
                return " {}['{}']".format(CONCENTRATIONS_DICT, parsed_id)
            else:
                raise ValueError("'{}' not a known specie in rate law".format(
                    parsed_id, rate_law_expression))

        if getattr(rate_law, 'equation', None) is None:
            return
        rate_law_equation = rate_law.equation
        if '__' in rate_law_equation.expression:
            raise ValueError("Security risk: rate law expression '{}' contains '__'.".format(
                rate_law_equation.expression))

        rate_law_expression = rate_law_equation.expression
        species_ids = set([specie.serialize() for specie in species])
        
        # Rate laws should be tokenized to properly construct a Python expression.
        # Simple pattern matching in which one specie name matches the suffix of another will fail.
        # Consider this string replace:
        #   py_expression = py_expression.replace(id, "concentrations['{}']".format(id))
        #   If these are two species: AB[c], CAB[c], then replace would produce Cconcentrations['AB[c]']
        #   which could not be evaluated.
        g = tokenize.generate_tokens(cStringIO(rate_law_equation.expression).readline)
        tokens = [(toknum, tokval) for toknum, tokval, _, _, _ in g]
        result = []
        idx = 0
        while idx < len(tokens):
            if possible_specie_id(tokens[idx:idx+4]):
                result.append(
                    (token.NAME, convert_specie_name(tokens[idx:], species_ids, rate_law_expression)))
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
            model (:obj:`Model`): a `Model`

        Raises:
            ValueError: If a rate law cannot be transcoded
        '''
        for submodel in model.get_submodels():
            for reaction in submodel.reactions:
                for rate_law in reaction.rate_laws:
                    rate_law.equation.transcoded = RateLawUtils.transcode(rate_law, 
                        submodel.get_species())

    @staticmethod
    def eval_rate_laws(reaction, concentrations):
        '''Evaluate a reaction's rate laws at the given species concentrations

        Args:
            reaction (:obj:`Reaction`): a Reaction instance
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
            transcoded_reaction = rate_law.equation.transcoded

            local_ns = {func.__name__: func for func in RateLawEquation.Meta.valid_functions}
            if hasattr(rate_law, 'k_cat') and not isnan(rate_law.k_cat):
                local_ns['k_cat'] = rate_law.k_cat
            if hasattr(rate_law, 'k_m') and not isnan(rate_law.k_m):
                local_ns['k_m'] = rate_law.k_m
            local_ns[CONCENTRATIONS_DICT] = concentrations

            try:
                rates.append(eval(transcoded_reaction, {}, local_ns))
            except SyntaxError as error:
                raise ValueError("Error: reaction '{}' has syntax error in transcoded rate law '{}'.".format(
                    reaction.id, transcoded_reaction))
            except NameError as error:
                raise NameError("Error: NameError in transcoded rate law '{}' of reaction '{}': '{}'".format(
                    transcoded_reaction, reaction.id, error))
            except Exception as error:
                raise Exception("Error: error in transcoded rate law '{}' of reaction '{}': '{}'".format(
                    transcoded_reaction, reaction.id, error))
        return rates
