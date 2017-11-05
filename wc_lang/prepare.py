""" Prepare a WC model for further processing, such as export or simulation

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-10-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

from math import ceil, floor, exp, log, log10, isnan
from warnings import warn
import ast

from obj_model import utils
from wc_utils.util.list import difference
from wc_lang.core import (SubmodelAlgorithm, Model, ObjectiveFunction, SpeciesType, SpeciesTypeType,
    Species, Concentration, Compartment, Reaction, ReactionParticipant, RateLawEquation, BiomassReaction)
from wc_lang.rate_law_utils import RateLawUtils

# configuration
from wc_utils.config.core import ConfigManager
from wc_lang.config import paths as config_paths_wc_lang
config_wc_lang = \
    ConfigManager(config_paths_wc_lang.core).get_config()['wc_lang']


class PrepareModel(object):
    '''Statically prepare a model

    `Models` which validate usually lack data needed to use them. `PrepareModel` automates
    the addition of default and statically computed data to a `Model`.

    Currently added data:
        Missing concentrations
        Fill gaps in dFBA submodel reaction networks
        Ensure that dFBA submodels have objective functions
        Apply default flux bounds to the reactions in dFBA submodels
    '''
    def __init__(self, model):
        self.model = model

    def run(self):
        """ Statically prepare a model by executing all `Prepare` methods.
        """
        for submodel in self.model.get_submodels():
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                reactions_created = self.fill_dfba_submodel_reaction_gaps(submodel)
                warn("{} gap filling reactions created for submodel '{}'.".format(reactions_created,
                    submodel.name))
                self.confirm_dfba_submodel_obj_func(submodel)
                (min_bounds_set, max_bounds_set) = self.apply_default_dfba_submodel_flux_bounds(submodel)
                warn("{} minimum and {} maximum default flux bounds set for submodel '{}'.".format(
                    min_bounds_set, max_bounds_set, submodel.name))
                try:
                    (reactions, biomass_reactions) = self.parse_dfba_submodel_obj_func(submodel)
                    PrepareModel.assign_linear_objective_fn(submodel, reactions, biomass_reactions)
                    submodel.objective_function.linear = True
                except Exception as e:
                    submodel.objective_function.linear = False
                    warn("Submodel '{}' has non-linear objective function '{}'.".format(submodel.name,
                        submodel.objective_function.expression))

        self.init_concentrations()

    @staticmethod
    def assign_linear_objective_fn(submodel, reactions, biomass_reactions):
        '''Assign a linear objective function to a submodel

        Assign a linear objective function parsed by `parse_dfba_submodel_obj_func` to a submodel's
        attributes.

        Args:
            submodel (`Submodel`): a dFBA submodel
            reactions (:obj:`list` of (`float`, `str`)): list of (coeff, id) pairs for reactions
            biomass_reactions (:obj:`list` of (`float`, `str`)): list of (coeff, id) pairs for
                biomass reactions

        Raises:
            ValueError: if `submodel` is not a dFBA submodel
        '''
        of = submodel.objective_function
        of.reactions = [Reaction.objects.get_one(id=id) for coeff,id in reactions]
        of.reaction_coefficients = [coeff for coeff,id in reactions]
        of.biomass_reactions = [BiomassReaction.objects.get_one(id=id) for coeff,id in biomass_reactions]
        of.biomass_reaction_coefficients = [coeff for coeff,id in biomass_reactions]

    def fill_dfba_submodel_reaction_gaps(self, submodel):
        '''Create reactions to fill gaps in a dFBA submodel's reaction network.

        Considering only reactions used by this submodel, generate gap filling reactions
        that consume species that are not consumed and produce species that are not produced. These
        reactions enable an FBA solver to obtain a steady state solution.

        Algorithm:
            S = the set of all species used by submodel
            generate a "-> s" reaction for each s in S that is not produced by a reaction in submodel
            generate an "s ->" reaction for each s in S that is not consumed by a reaction in submodel

        Args:
            submodel (`Submodel`): a DFBA submodel

        Raises:
            ValueError: if `submodel` is not a dFBA submodel

        Returns:
            (:obj:`int`): the number of reactions created
        '''
        if submodel.algorithm != SubmodelAlgorithm.dfba:
            raise ValueError("submodel '{}' not a dfba submodel".format(submodel.name))

        reaction_number = 1

        species = submodel.get_species()
        species_not_produced = set(species)
        species_not_consumed = set(species)
        for rxn in submodel.reactions:
            if rxn.reversible:
                for part in rxn.participants:
                    species_not_produced.discard(part.species)
                    species_not_consumed.discard(part.species)
            else:
                for part in rxn.participants:
                    if part.coefficient<0:
                        species_not_consumed.discard(part.species)
                    elif 0<part.coefficient:
                        species_not_produced.discard(part.species)

        if species_not_produced & species_not_consumed:
            warn("some species in submodel '{}' are neither produced nor consumed: {}".format(
                submodel.id,
                sorted([s.id() for s in species_not_produced & species_not_consumed])))

        GAP_FILLING_RXN_ID_PREFIX = config_wc_lang['GAP_FILLING_RXN_ID_PREFIX']
        GAP_FILLING_RXN_NAME_PREFIX = config_wc_lang['GAP_FILLING_RXN_NAME_PREFIX']
        for specie in species_not_produced:
            # generate a "-> specie" reaction
            new_rxn = submodel.reactions.create(
                id = "{}_{}".format(GAP_FILLING_RXN_ID_PREFIX, reaction_number),
                name = "{}_{}".format(GAP_FILLING_RXN_NAME_PREFIX, reaction_number))
            reaction_number += 1
            new_rxn.participants.create(species=specie, coefficient = 1)

        for specie in species_not_consumed:
            # generate a "specie ->" reaction
            new_rxn = submodel.reactions.create(
                id = "{}_{}".format(GAP_FILLING_RXN_ID_PREFIX, reaction_number),
                name = "{}_{}".format(GAP_FILLING_RXN_NAME_PREFIX, reaction_number))
            reaction_number += 1
            new_rxn.participants.create(species=specie, coefficient = -1)

        return reaction_number-1

    def confirm_dfba_submodel_obj_func(self, submodel):
        '''Ensure that a dFBA submodel has an objective function

        If the submodel definition does not provide an objective function, then use the
        biomass reaction.

        Args:
            submodel (`Submodel`): a dFBA submodel

        Raises:
            ValueError: if `submodel` is not a dFBA submodel
            ValueError: if `submodel` cannot use its biomass reaction '{}' as an objective function
        '''
        if submodel.algorithm != SubmodelAlgorithm.dfba:
            raise ValueError("submodel '{}' not a dfba submodel".format(submodel.name))

        if not submodel.objective_function is None:
            submodel.objective_function.reaction_coefficients = []
            submodel.objective_function.biomass_reaction_coefficients = []
            return None

        # use the biomass reaction as the objective, because no objective function is specified
        of = ObjectiveFunction(expression=submodel.biomass_reaction.id,
            reactions=[],
            biomass_reactions=[submodel.biomass_reaction])

        of.reaction_coefficients = []
        of.biomass_reaction_coefficients = [1.0]
        submodel.objective_function = of
        return None

    def parse_dfba_submodel_obj_func(self, submodel):
        '''Parse a dFBA submodel's objective function into a linear function of reaction fluxes

        The SBML FBC only handles objectives that are a linear function of reaction fluxes. This method
        uses Python's parser to parse an objective function.

        The general form for an objective is :math:`c_1*id_1 + c_2*id_2 + ... + c_n*id_n`,
        where :math:`c_i` is a numerical coefficient and :math:`id_i` is a reaction id. The ids may
        represent reactions or biomass reactions.
        Coefficients may also appear after an id, as in :math:`id_j*c_j`. Coefficients equal to 1.0
        are not needed. And negative coefficients are supported.

        Args:
            submodel (`Submodel`): a dFBA submodel

        Returns:
            :obj:`(`list`, `list`)`: a pair of lists representing the objective's linear form;
                (`list` of (coeff, id) pairs for reactions, `list` of (coeff, id) pairs for biomass
                reactions)

        Raises:
            ValueError: if `submodel` is not a dFBA submodel
            ValueError: if `submodel.objective_function` is not a legal python expression, does not
                have the form above, is not a linear function of reaction ids, uses an unknown
                reaction id, or uses an id multiple times
        '''
        if submodel.algorithm != SubmodelAlgorithm.dfba:
            raise ValueError("submodel '{}' not a dfba submodel".format(submodel.name))

        def proc_mult(node, linear_expr):
            ''' Process a Mult node in the ast.

            Append the Mult node's coefficient and reaction id to `linear_expr`.

            Args:
                node (:obj:`ast.BinOp`): an ast binary operation that uses multiplication
                linear_expr (:obj:`list` of `tuple`): pairs of (coefficient, reaction_id)

            Raises:
                :obj:`ValueError`: if the Mult node does not have one Name and one Num (which may be negative)
            '''
            nums = []
            names = []
            sign = 1.0
            for element in [node.left, node.right]:
                if isinstance(element, ast.Num):
                    nums.append(element)
                if isinstance(element, ast.UnaryOp) and isinstance(element.op, ast.USub):
                    # the coefficient is negative
                    sign = -1.0
                    nums.append(element.operand)
                if isinstance(element, ast.Name):
                    names.append(element)
            if not (len(nums)==1 and len(names)==1):
                raise ValueError("bad Mult")
            linear_expr.append((sign*nums[0].n, names[0].id))

        def proc_add(node, linear_expr):
            ''' Process an Add node in the ast.

            Append the Add node's coefficient(s) and reaction id(s) to `linear_expr`.

            Args:
                node (:obj:`ast.BinOp`): an ast binary operation that uses addition
                linear_expr (:obj:`list` of `tuple`): pairs of (coefficient, reaction_id)

            Raises:
                :obj:`ValueError`: if the Add node does not have a total of 2 Names, Mults, and Adds.
            '''
            names = []
            mults = []
            adds = 0
            for element in [node.left, node.right]:
                if isinstance(element, ast.Name):
                    names.append(element)
                if isinstance(element, ast.BinOp):
                    if isinstance(element.op, ast.Mult):
                        mults.append(element)
                    if isinstance(element.op, ast.Add):
                        adds += 1
            if len(names) + len(mults) + adds != 2:
                raise ValueError("bad Add")
            # A Name that's not in a mult. op. is multiplied by 1.0 by default
            # An Add may contain 2 of them
            for name in names:
                linear_expr.append((1.0, name.id))

        linear_expr = []    # list of (coeff, reaction_id)
        objective_function = submodel.objective_function
        objective_function.expression = objective_function.expression.strip()
        expected_nodes = (ast.Add, ast.Expression, ast.Load, ast.Mult, ast.Num, ast.USub,
            ast.UnaryOp, ast.Name)
        try:
            for node in ast.walk(ast.parse(objective_function.expression, mode='eval')):
                try:
                    # if linear_expr is empty then an ast.Name is the entire expression
                    if isinstance(node, ast.Name) and not linear_expr:
                        linear_expr.append((1.0, node.id))
                    elif isinstance(node, ast.BinOp):
                        if isinstance(node.op, ast.Mult):
                            proc_mult(node, linear_expr)
                        elif isinstance(node.op, ast.Add):
                            proc_add(node, linear_expr)
                    elif isinstance(node, expected_nodes):
                        continue
                    else:
                        raise ValueError()
                except ValueError:
                    raise ValueError("Cannot parse objective function '{}' as a linear function of "
                        "reaction ids.".format(objective_function.expression))
        except Exception as e:
            raise ValueError("Cannot parse objective function '{}'.".format(objective_function.expression))

        # error if multiple uses of a reaction in an objective function
        seen = set()
        dupes = []
        for id in [id for coeff,id in linear_expr]:
            if id in seen and id not in dupes:
                dupes.append(id)
            seen.add(id)
        if dupes:
            raise ValueError("Multiple uses of '{}' in objective function '{}'.".format(dupes,
                objective_function.expression))
        reactions = []
        biomass_reactions = []

        for coeff,id in linear_expr:

            if Reaction.objects.get_one(id=id):
                reactions.append((coeff,id),)
                continue

            if BiomassReaction.objects.get_one(id=id):
                biomass_reactions.append((coeff,id),)
                continue

            raise ValueError("Unknown reaction id '{}' in objective function '{}'.".format(id,
                objective_function.expression))
        return (reactions, biomass_reactions)

    def apply_default_dfba_submodel_flux_bounds(self, submodel):
        ''' Apply default flux bounds to a dFBA submodel's reactions

        The FBA optimizer needs min and max flux bounds for each dFBA submodel reaction.
        If some reactions lack bounds and default bounds are provided in a config file,
        then apply the defaults to the reactions.
        Specifically, min and max default bounds are applied as follows:

            reversible reactions:

              * min_flux = -default_max_flux_bound
              * max_flux = default_max_flux_bound

            irreversible reactions:

              * min_flux = default_min_flux_bound
              * max_flux = default_max_flux_bound

        Args:
            submodel (`Submodel`): a dFBA submodel

        Raises:
            ValueError: if `submodel` is not a dFBA submodel

        Returns:
            :obj:`tuple` of (`int`, `int`): counts of min and max flux bounds set to the default
        '''
        if submodel.algorithm != SubmodelAlgorithm.dfba:
            raise ValueError("submodel '{}' not a dfba submodel".format(submodel.name))

        need_default_flux_bounds = False
        for rxn in submodel.reactions:
            need_default_flux_bounds = need_default_flux_bounds or isnan(rxn.min_flux) or isnan(rxn.max_flux)
        if not need_default_flux_bounds:
            # all reactions have flux bounds
            return (0,0)

        # Are default flux bounds available? They cannot be negative.
        try:
            default_min_flux_bound = config_wc_lang['default_min_flux_bound']
            default_max_flux_bound = config_wc_lang['default_max_flux_bound']
        except KeyError as e:
            raise ValueError("cannot obtain default_min_flux_bound and default_max_flux_bound=")
        if not 0 <= default_min_flux_bound <= default_max_flux_bound:
            raise ValueError("default flux bounds violate 0 <= default_min_flux_bound <= default_max_flux_bound:\n"
            "default_min_flux_bound={}; default_max_flux_bound={}".format(default_min_flux_bound,
                default_max_flux_bound))

        # Apply default flux bounds to reactions in submodel
        num_default_min_flux_bounds = 0
        num_default_max_flux_bounds = 0
        for rxn in submodel.reactions:
            if isnan(rxn.min_flux):
                num_default_min_flux_bounds += 1
                if rxn.reversible:
                    rxn.min_flux = -default_max_flux_bound
                else:
                    rxn.min_flux = default_min_flux_bound
            if isnan(rxn.max_flux):
                num_default_max_flux_bounds += 1
                rxn.max_flux = default_max_flux_bound
        return (num_default_min_flux_bounds, num_default_max_flux_bounds)

    def init_concentrations(self):
        """ Initialize missing concentration values to 0 """
        for specie in self.model.get_species():
            if specie.concentration is None:
                warn("setting concentration for {} to 0.0".format(specie.id()))
                specie.concentrations = Concentration(species=specie, value=0.0)

# TODO: fix doc string formatting
class CheckModel(object):
    '''Statically check a model

    A `Model` which validates may fail to satisfy global properties that must hold for the `Model`
    to be used. `CheckModel` evaluates these properties.

    Currently checked properties:
        DFBA submodels contain a biomass reaction and an objective function
        Rate laws transcode and evaluate without error
        All reactants in each submodel's reactions are in the submodel's compartment

    Other properties to check:
        The model does not contain dead-end species which are only consumed or produced
        Reactions are balanced
        Reactions in dynamic submodels contain fully specified rate laws
        All Species used in reactions have concentration values
        Consider the reactions modeled by a submodel -- all modifier species used by the rate laws
            for the reactions participate in at least one reaction in the submodel
        Ensure that Reaction and BiomassReaction ids don't overlap; can then simplify ObjectiveFunction.deserialize()
        # TODO: implement these, and expand the list of properties
    '''
    def __init__(self, model):
        self.model = model

    def run(self):
        self.errors = []
        for submodel in self.model.get_submodels():
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                self.errors.extend(self.check_dfba_submodel(submodel))
            if submodel.algorithm in [SubmodelAlgorithm.ssa, SubmodelAlgorithm.ode]:
                self.errors.extend(self.check_dynamic_submodel(submodel))
        self.errors.extend(self.transcode_and_check_rate_law_equations())
        if self.errors:
            raise ValueError('\n'.join(self.errors))

    def check_dfba_submodel(self, submodel):
        '''Check the inputs to a DFBA submodel

        Ensure that:
            * All regular DFBA reactions have min flux and max flux with appropriate values
            * The DFBA submodel contains a biomass reaction and an objective function

        Args:
            submodel (`Submodel`): a DFBA submodel

        Returns:
            :obj:`list` of `str`: if no errors, returns an empty `list`, otherwise a `list` of
                error messages
        '''
        errors = []
        for reaction in submodel.reactions:
            for attr in ['min_flux', 'max_flux']:
                if not hasattr(reaction, attr) or isnan(getattr(reaction, attr)):
                    errors.append("Error: no {} for reaction '{}' in submodel '{}'".format(
                        attr, reaction.name, submodel.name))
                    continue

            if hasattr(reaction, 'min_flux') and hasattr(reaction, 'max_flux'):
                if reaction.max_flux < reaction.min_flux:
                    errors.append("Error: max_flux < min_flux ({} < {}) for reaction '{}' in submodel '{}'".format(
                        reaction.max_flux, reaction.min_flux, reaction.name, submodel.name))
                if reaction.reversible and 0 < reaction.min_flux:
                    errors.append("Error: 0 < min_flux ({}) for reversible reaction '{}' in submodel '{}'".format(
                        reaction.min_flux, reaction.name, submodel.name))

        if submodel.biomass_reaction is None or not submodel.biomass_reaction.biomass_components:
            errors.append("Error: submodel '{}' uses dfba but lacks a biomass reaction".format(submodel.name))

        if submodel.objective_function is None:
            errors.append("Error: submodel '{}' uses dfba but lacks an objective function".format(submodel.name))

        return errors

    def check_dynamic_submodel(self, submodel):
        '''Check the inputs to a dynamic submodel

        Ensure that:
            * All reactions have rate laws for the appropriate directions

        Args:
            submodel (`Submodel`): a dynamic (SSA or ODE) submodel

        Returns:
            :obj:`list` of :obj:`str` if no errors, returns an empty `list`, otherwise a `list` of
                error messages
        '''
        errors = []
        for reaction in submodel.reactions:
            direction_types = set()
            for rate_law in reaction.rate_laws:
                direction_types.add(rate_law.direction.name)
            if not direction_types:
                errors.append("Error: reaction '{}' in submodel '{}' has no "
                    "rate law specified".format(reaction.name, submodel.name))
            if reaction.reversible:     # reversible is redundant with a reaction's rate laws
                if direction_types.symmetric_difference(set(('forward', 'backward'))):
                    errors.append("Error: reaction '{}' in submodel '{}' is reversible but has only "
                        "a '{}' rate law specified".format(reaction.name, submodel.name,
                        direction_types.pop()))
            else:
                if direction_types.symmetric_difference(set(('forward',))):
                    errors.append("Error: reaction '{}' in submodel '{}' is not reversible but has "
                        "a 'backward' rate law specified".format(reaction.name, submodel.name))
        return errors

    def transcode_and_check_rate_law_equations(self):
        '''Transcode and evaluate all rate law equations in a model

        Ensure that all rate law equations can be transcoded and evaluated.

        Returns:
            :obj:`list` of `str`: if no errors, returns an empty `list`, otherwise a `list` of
                error messages
        '''
        errors = []
        species = self.model.get_species()
        species_concentrations = {}
        for concentration in self.model.get_concentrations():
            species_concentrations[concentration.species.serialize()] = concentration.value
        for reaction in self.model.get_reactions():
            for rate_law in reaction.rate_laws:
                if getattr(rate_law, 'equation', None) is None:
                    continue
                try:
                    rate_law.equation.transcoded = RateLawUtils.transcode(rate_law.equation, species)
                except Exception as error:
                    errors.append(str(error))
            try:
                rates = RateLawUtils.eval_reaction_rate_laws(reaction, species_concentrations)
            except Exception as error:
                errors.append(str(error))
        return errors

    def verify_reactant_compartments(self):
        '''Verify that all reactants in each submodel's reactions are in the submodel's compartment

        Returns:
            :obj:`list` of `str`: if no errors, returns an empty `list`, otherwise a `list` of
                error messages
        '''
        errors = []
        for lang_submodel in self.model.get_submodels():
            compartment = lang_submodel.compartment
            if compartment is None:
                errors.append("submodel '{}' must contain a compartment attribute".format(lang_submodel.id))
                continue
            for reaction in lang_submodel.reactions:
                for participant in reaction.participants:
                    if participant.coefficient < 0:     # select reactants
                        if participant.species.compartment != compartment:
                            error = "submodel '{}' models compartment {}, but its reaction {} uses "\
                            "specie {} in another compartment: {}".format(lang_submodel.id,
                                compartment.id, reaction.id, participant.species.species_type.id,
                                participant.species.compartment.id)
                            errors.append(error)
        return errors
