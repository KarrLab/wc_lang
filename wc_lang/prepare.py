""" Prepare a WC model for further processing, such as export or simulation.

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-10-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

from math import ceil, floor, exp, log, log10, isnan
from warnings import warn
import ast
import networkx as nx

from obj_model import utils
from wc_utils.util.list import difference
from obj_model.utils import get_component_by_id
from wc_lang import (SubmodelAlgorithm, Model, ObjectiveFunction, SpeciesType, SpeciesTypeType,
                     Species, Concentration, Compartment, Reaction, SpeciesCoefficient, RateLawEquation,
                     BiomassReaction, ObservableExpression, Observable, Function, FunctionExpression)
from wc_lang.expression_utils import RateLawUtils

# configuration
import wc_lang.config.core
config_wc_lang = wc_lang.config.core.get_config()['wc_lang']

EXTRACELLULAR_COMPARTMENT_ID = config_wc_lang['EXTRACELLULAR_COMPARTMENT_ID']

# TODO: distinguish between preparing and analyzing a model: analyze includes gap finding


class AnalyzeModel(object):
    """ Statically analyze a model

    `AnalyzeModel` performs static analysis of WC-lang models which are useful for constructing
    models.

    Current analyses:

        * Identify dead end species and reaction network gaps in dFBA submodels
    """

    def __init__(self, model):
        self.model = model

    def identify_dfba_submodel_rxn_gaps(self, submodel):
        """ Identify gaps in a dFBA submodel's reaction network

        Species that are not consumed or not produced indicate gaps in the reaction network.
        These can be found by a static analysis of the model. Reactions that use species that
        are not produced or produce species that are not consumed must eventually have zero flux.
        A reaction network can be reduced to a minimal network of reactions that can all
        have positive fluxes.

        Algorithm::

            all_gap_species = get_gap_species([])
            delta_gap_species = all_gap_species
            while delta_gap_species:
                all_gap_reactions = get_gap_reactions(all_gap_species)
                tmp_gap_species = all_gap_species
                all_gap_species = get_gap_species(all_gap_reactions)
                delta_gap_species = all_gap_species - tmp_gap_species
            return (all_gap_species, all_gap_reactions)

        Args:
            submodel (`Submodel`): a DFBA submodel

        Raises:
            ValueError: if `submodel` is not a dFBA submodel

        Returns:
            :obj:`tuple`:

                * :obj:`set` of :obj:`Species`: `Species` not in the minimal reaction network
                * :obj:`set` of :obj:`Reaction`: `Reaction`s not in the minimal reaction network
        """
        if submodel.algorithm != SubmodelAlgorithm.dfba:
            raise ValueError("submodel '{}' not a dfba submodel".format(submodel.name))

        all_dead_end_species = AnalyzeModel.find_dead_end_species(submodel, set())
        delta_dead_end_species = all_dead_end_species
        inactive_reactions = set()
        while any(delta_dead_end_species):
            inactive_reactions = AnalyzeModel.get_inactive_reactions(submodel, all_dead_end_species)
            tmp_not_consumed, tmp_not_produced = all_dead_end_species
            all_dead_end_species = AnalyzeModel.find_dead_end_species(submodel, inactive_reactions)
            all_not_consumed, all_not_produced = all_dead_end_species
            delta_dead_end_species = (all_not_consumed-tmp_not_consumed, all_not_produced-tmp_not_produced)
        return(all_dead_end_species, inactive_reactions)

    @staticmethod
    def find_dead_end_species(submodel, inactive_reactions):
        """ Find the dead end species in a reaction network

        Given a set of inactive reactions in submodel, determine species that are not consumed by
        any reaction, or are not produced by any reaction. Costs :math:`O(n*p)`, where :math:`n` is
        the number of reactions in `submodel` and :math:`p` is the maximum number of participants in
        a reaction.

        Args:
            submodel (`Submodel`): a DFBA submodel
            inactive_reactions (`set` of `Reaction`): the inactive reactions in `submodel`

        Returns:
            :obj:`tuple`:

                * :obj:`set` of :obj:`Species`: the species that are not consumed
                * :obj:`set` of :obj:`Species`: the species that are not produced
        """
        species = submodel.get_species()
        species_not_consumed = set(species)
        species_not_produced = set(species)
        for rxn in submodel.reactions:
            if rxn in inactive_reactions:
                continue
            if rxn.reversible:
                for part in rxn.participants:
                    species_not_consumed.discard(part.species)
                    species_not_produced.discard(part.species)
            else:
                for part in rxn.participants:
                    if part.coefficient < 0:
                        species_not_consumed.discard(part.species)
                    elif 0 < part.coefficient:
                        species_not_produced.discard(part.species)
        return(species_not_consumed, species_not_produced)

    @staticmethod
    def get_inactive_reactions(submodel, dead_end_species):
        """ Find the inactive reactions in a reaction network

        Given the dead end species in a reaction network, find the reactions that must eventually
        become inactive. Reactions that consume species which are not produced must become inactive.
        And reactions that produce species which are not consumed must become inactive to prevent
        the copy numbers of those species from growing without bound.
        Costs :math:`O(n*p)`, where :math:`n` is the number of reactions in `submodel` and :math:`p`
        is the maximum number of participants in a reaction.

        Args:
            submodel (:obj:`Submodel`): a DFBA submodel
            dead_end_species (:obj:`tuple`):

                * :obj:`set` of :obj:`Species`: the `Species` that are not consumed by any `Reaction` in `submodel`
                * :obj:`set` of :obj:`Species`: the `Species` that are not produced by any `Reaction` in `submodel`

        Returns:
            :obj:`set` of :obj:`Reaction`: the inactive reactions in `submodel`'s reaction network
        """
        species_not_consumed, species_not_produced = dead_end_species
        inactive_reactions = []
        for rxn in submodel.reactions:
            for part in rxn.participants:
                if (part.species in species_not_consumed or
                        part.species in species_not_produced):
                    inactive_reactions.append(rxn)
                    break
        return inactive_reactions

    @staticmethod
    def digraph_of_rxn_network(submodel):
        """ Create a NetworkX network representing the reaction network in `submodel`

        To leverage the algorithms in NetworkX, map a reaction network onto a NetworkX
        directed graph.
        The digraph is bipartite, with `Reaction` and `Species` nodes. A reaction is represented
        a Reaction node, with an edge from each reactant Species node to the Reaction node, and
        an edge from the Reaction node to each product Species node.

        Args:
            submodel (:obj:`Submodel`): a DFBA submodel

        Returns:
            :obj:`DiGraph`: a NetworkX directed graph representing `submodel`'s reaction network
        """
        digraph = nx.DiGraph()

        # make network of obj_model.Model instances
        for specie in submodel.get_species():
            digraph.add_node(specie)
        for rxn in submodel.reactions:
            digraph.add_node(rxn)
            for participant in rxn.participants:
                part = participant.species
                if participant.coefficient < 0:
                    # reactant
                    digraph.add_edge(part, rxn)
                elif 0 < participant.coefficient:
                    # product
                    digraph.add_edge(rxn, part)
            if rxn.reversible:
                for participant in rxn.participants:
                    part = participant.species
                    if participant.coefficient < 0:
                        # product
                        digraph.add_edge(rxn, part)
                    elif 0 < participant.coefficient:
                        # reactant
                        digraph.add_edge(part, rxn)
        return digraph

    @staticmethod
    def path_bounds_analysis(submodel):
        """ Perform path bounds analysis on `submodel`

        To be adequately constrained, a dFBA metabolic model should have the property that each path
        from an extracellular species to a component in the objective function contains at least
        one reaction constrained by a finite flux upper bound.

        Analyze the reaction network in `submodel` and return all paths from extracellular species
        to objective function components that lack a finite flux upper bound.

        Args:
            submodel (:obj:`Submodel`): a DFBA submodel

        Returns:
            :obj:`dict` of `list` of `list` of :obj:: paths from extracellular species to objective
            function components that lack a finite flux upper bound. Keys in the `dict` are the ids
            of extracellular species; the corresponding values contain the unbounded paths for the
            extracellular species, as returned by `unbounded_paths`.
        """
        # todo: symmetrically, report reactions not on any path from ex species to obj fun components
        digraph = AnalyzeModel.digraph_of_rxn_network(submodel)
        obj_fn_species = submodel.objective_function.get_products()
        ex_compartment = submodel.model.compartments.get_one(id=EXTRACELLULAR_COMPARTMENT_ID)
        ex_species = submodel.get_species(compartment=ex_compartment)
        all_unbounded_paths = dict()
        for ex_specie in ex_species:
            paths = AnalyzeModel.unbounded_paths(digraph, ex_specie, obj_fn_species)
            all_unbounded_paths[ex_specie.id()] = paths
        return all_unbounded_paths

    # todo: replace the constant in min_non_finite_ub=1000.0
    @staticmethod
    def unbounded_paths(rxn_network, ex_species, obj_fn_species, min_non_finite_ub=1000.0):
        """ Find the unbounded paths from an extracellular species to some objective function species

        Return all paths in a reaction network that lack a finite flux upper bound
        and go from `ex_species` to an objective function component.

        Args:
            rxn_network (:obj:`DiGraph`): a NetworkX directed graph representing a reaction network,
            created by `digraph_of_rxn_network`
            ex_species (:obj:`Species`): an extracellular `Species` that is a node in `rxn_network`
            obj_fn_species (:obj:`list` of :obj:`Species`): objective function `Species` that are
            also nodes in `rxn_network`
            finite_upper_bound_limit (:obj:`float`, optional): the maximum value of a finite flux
            upper bound
            min_non_finite_ub (:obj:`float`, optional): flux upper bounds less than `min_non_finite_ub`
            are considered finite

        Returns:
            :obj:`list` of `list` of :obj:: a list of the reaction paths from `ex_species`
            to objective function components that lack a finite flux upper bound.
            A path is a list of `Species`, `Reaction`, `Species`, ..., `Species`, starting with
            `ex_species` and ending with an objective function component.
        """
        unbounded_paths = list()
        if not isinstance(ex_species, Species):
            raise ValueError("'ex_species' should be a Species instance, but it is a {}".format(
                type(ex_species).__name__))
        for of_specie in obj_fn_species:
            if not isinstance(of_specie, Species):
                raise ValueError("elements of 'obj_fn_species' should be Species instances, but one is a {}".format(
                    type(of_specie).__name__))
            for path in nx.all_simple_paths(rxn_network, source=ex_species, target=of_specie):
                # path is a list of Species, Reaction, ..., Species
                bounded = False
                for i in range(1, len(path), 2):
                    rxn = path[i]
                    if rxn.max_flux < min_non_finite_ub:
                        bounded = True
                        break
                if not bounded:
                    unbounded_paths.append(path)
        return unbounded_paths


class PrepareModel(object):
    """ Statically prepare a model

    `Models` which validate usually lack data needed to use them. `PrepareModel` automates
    the addition of default and statically computed data to a `Model`.

    Currently added data:

        * Missing concentrations
        * Create implicit exchange reactions for dFBA submodels
        * Ensure that dFBA submodels have objective functions
        * Apply default flux bounds to the reactions in dFBA submodels
    """

    def __init__(self, model):
        self.model = model

    def run(self):
        """ Statically prepare a model by executing all `Prepare` methods.
        """
        for submodel in self.model.get_submodels():
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                reactions_created = self.create_dfba_exchange_rxns(submodel,
                                                                   EXTRACELLULAR_COMPARTMENT_ID)
                warn("{} exchange reactions created for submodel '{}'.".format(reactions_created,
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

    def create_dfba_exchange_rxns(self, submodel, extracellular_compartment_id):
        """ Create exchange reactions for a dFBA submodel's reaction network.

        To represent FBA's mathematical assumption that it models a closed system, create
        'implicit' forward exchange reactions that synthesize all extracellular metabolites.

        # TODO: To model how other pathways consume metabolites generated by metabolism, create 'implicit'
        reactions which exchange these metabolites between a dFBA metabolism submodel and the other
        pathway(s)/submodel(s).

        Algorithm to synthesize extracellular metabolites::

            E = the set of all extracellular metabolites used by the submodel
            generate a "-> e" reaction for each e in E in the submodel

        Args:
            submodel (`Submodel`): a DFBA submodel
            extracellular_compartment_id (`str`): the id of the extracellular compartment

        Raises:
            ValueError: if `submodel` is not a dFBA submodel

        Returns:
            :obj:`int`: the number of reactions created
        """
        if submodel.algorithm != SubmodelAlgorithm.dfba:
            raise ValueError("submodel '{}' not a dfba submodel".format(submodel.name))

        reaction_number = 1

        for specie in submodel.get_species():
            if specie.compartment.id == extracellular_compartment_id:

                EXCHANGE_RXN_ID_PREFIX = config_wc_lang['EXCHANGE_RXN_ID_PREFIX']
                EXCHANGE_RXN_NAME_PREFIX = config_wc_lang['EXCHANGE_RXN_NAME_PREFIX']
                # generate a "-> specie" reaction
                new_rxn = submodel.reactions.create(
                    id="{}_{}".format(EXCHANGE_RXN_ID_PREFIX, reaction_number),
                    name="{}_{}".format(EXCHANGE_RXN_NAME_PREFIX, reaction_number),
                    reversible=False,
                    min_flux=-float('inf'),
                    max_flux=float('inf'))
                reaction_number += 1
                new_rxn.participants.create(species=specie, coefficient=1)

        return reaction_number-1

    def confirm_dfba_submodel_obj_func(self, submodel):
        """ Ensure that a dFBA submodel has an objective function

        If the submodel definition does not provide an objective function, then use the
        biomass reaction.

        Args:
            submodel (`Submodel`): a dFBA submodel

        Raises:
            ValueError: if `submodel` is not a dFBA submodel
            ValueError: if `submodel` cannot use its biomass reaction '{}' as an objective function
        """
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
        """ Parse a dFBA submodel's objective function into a linear function of reaction fluxes

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
            :obj:`tuple`: a pair of lists representing the objective's linear form

                * obj:`list`: (coeff, id) pairs for reactions
                * obj:`list`: (coeff, id) pairs for biomass reactions

        Raises:
            ValueError: if `submodel` is not a dFBA submodel
            ValueError: if `submodel.objective_function` is not a legal python expression, does not
                have the form above, is not a linear function of reaction ids, uses an unknown
                reaction id, or uses an id multiple times
        """
        if submodel.algorithm != SubmodelAlgorithm.dfba:
            raise ValueError("submodel '{}' not a dfba submodel".format(submodel.name))

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
                            self._proc_mult(node, linear_expr)
                        elif isinstance(node.op, ast.Add):
                            self._proc_add(node, linear_expr)
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
        for id in [id for coeff, id in linear_expr]:
            if id in seen and id not in dupes:
                dupes.append(id)
            seen.add(id)
        if dupes:
            raise ValueError("Multiple uses of '{}' in objective function '{}'.".format(dupes,
                                                                                        objective_function.expression))
        reactions = []
        biomass_reactions = []

        for coeff, id in linear_expr:

            reaction = get_component_by_id(submodel.model.get_reactions(), id)
            if reaction:
                reactions.append((coeff, id),)
                continue

            biomass_reaction = get_component_by_id(submodel.model.get_biomass_reactions(), id)
            if biomass_reaction:
                biomass_reactions.append((coeff, id),)
                continue

            raise ValueError("Unknown reaction or biomass reaction id '{}' in objective function '{}'.".format(
                id, objective_function.expression))
        return (reactions, biomass_reactions)

    @staticmethod
    def _proc_mult(node, linear_expr):
        """ Process a Mult node in the ast.

        Append the Mult node's coefficient and reaction id to `linear_expr`.

        Args:
            node (:obj:`ast.BinOp`): an ast binary operation that uses multiplication
            linear_expr (:obj:`list` of `tuple`): pairs of (coefficient, reaction_id)

        Raises:
            :obj:`ValueError`: if the Mult node does not have one Name and one Num (which may be negative)
        """
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
        if not (len(nums) == 1 and len(names) == 1):
            raise ValueError("bad Mult")
        linear_expr.append((sign*nums[0].n, names[0].id))

    @staticmethod
    def _proc_add(node, linear_expr):
        """ Process an Add node in the ast.

        Append the Add node's coefficient(s) and reaction id(s) to `linear_expr`.

        Args:
            node (:obj:`ast.BinOp`): an ast binary operation that uses addition
            linear_expr (:obj:`list` of `tuple`): pairs of (coefficient, reaction_id)

        Raises:
            :obj:`ValueError`: if the Add node does not have a total of 2 Names, Mults, and Adds.
        """
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

    @staticmethod
    def assign_linear_objective_fn(submodel, reactions, biomass_reactions):
        """ Assign a linear objective function to a submodel

        Assign a linear objective function parsed by `parse_dfba_submodel_obj_func` to a submodel's
        attributes.

        Args:
            submodel (`Submodel`): a dFBA submodel
            reactions (:obj:`list` of (`float`, `str`)): list of (coeff, id) pairs for reactions
            biomass_reactions (:obj:`list` of (`float`, `str`)): list of (coeff, id) pairs for
                biomass reactions

        Raises:
            ValueError: if `submodel` is not a dFBA submodel
        """
        of = submodel.objective_function
        of.reactions = [get_component_by_id(submodel.model.get_reactions(), id) for coeff, id in reactions]
        of.reaction_coefficients = [coeff for coeff, id in reactions]
        of.biomass_reactions = [get_component_by_id(submodel.model.get_biomass_reactions(), id)
                                for coeff, id in biomass_reactions]
        of.biomass_reaction_coefficients = [coeff for coeff, id in biomass_reactions]

    def apply_default_dfba_submodel_flux_bounds(self, submodel):
        """ Apply default flux bounds to a dFBA submodel's reactions

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
            submodel (:obj:`Submodel`): a dFBA submodel

        Raises:
            ValueError: if `submodel` is not a dFBA submodel

        Returns:
            :obj:`tuple`:

                * obj:`int`: number of min flux bounds set to the default
                * obj:`int`: number of max flux bounds set to the default
        """
        if submodel.algorithm != SubmodelAlgorithm.dfba:
            raise ValueError("submodel '{}' not a dfba submodel".format(submodel.name))

        need_default_flux_bounds = False
        for rxn in submodel.reactions:
            need_default_flux_bounds = need_default_flux_bounds or isnan(rxn.min_flux) or isnan(rxn.max_flux)
        if not need_default_flux_bounds:
            # all reactions have flux bounds
            return (0, 0)

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


class CheckModel(object):
    """ Statically check a model

    A `Model` which validates in `wc_lang` may fail to satisfy global properties that must hold for
    the `Model` to be used. `CheckModel` evaluates these properties.

    Currently checked properties:

        * DFBA submodels contain a biomass reaction and an objective function
        * Rate laws transcode and evaluate without error
        * All reactants in each submodel's reactions are in the submodel's compartment
        * All species types have positive molecular weights
        * The network of `Observable` depencencies is acyclic

    Other properties to be checked:

        * The model does not contain dead-end species which are only consumed or produced
        * The model uses water consistently - either in all compartments or in none
        * Reactions are balanced
        * Reactions in dynamic submodels contain fully specified rate laws
        * A reaction's rate laws uses only species that are in the reaction's reactants
        * All Species used in reactions have concentration values
        * Consider the reactions modeled by a submodel -- all modifier species used by the rate laws
          for the reactions participate in at least one reaction in the submodel
        * Ensure that Reaction and BiomassReaction ids don't overlap; can then simplify ObjectiveFunction.deserialize()

    # TODO: implement these, and expand the list of properties
    """

    def __init__(self, model):
        self.model = model

    def run(self):
        """ Run all tests in `CheckModel`

        Raises:
            :obj:`ValueError`: if any of the tests return an error
        """
        self.errors = []
        for submodel in self.model.get_submodels():
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                self.errors.extend(self.check_dfba_submodel(submodel))
            if submodel.algorithm in [SubmodelAlgorithm.ssa, SubmodelAlgorithm.ode]:
                self.errors.extend(self.check_dynamic_submodel(submodel))
        self.errors.extend(self.transcode_and_check_rate_law_equations())
        self.errors.extend(self.verify_species_types())
        self.errors.extend(self.verify_acyclic_dependencies([Observable]))
        if self.errors:
            raise ValueError('\n'.join(self.errors))

    def check_dfba_submodel(self, submodel):
        """ Check the inputs to a DFBA submodel

        Ensure that:

            * All DFBA reactions have min flux and max flux with appropriate values
            * The DFBA submodel contains a biomass reaction and an objective function
            * All species used in biomass reactions are defined

        Args:
            submodel (`Submodel`): a DFBA submodel

        Returns:
            :obj:`list` of :obj:`str`: if no errors, returns an empty `list`; otherwise a `list` of
            error messages
        """
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

        if submodel.objective_function is None:
            errors.append("Error: submodel '{}' uses dfba but lacks an objective function".format(submodel.name))

        if submodel.biomass_reaction is None or not submodel.biomass_reaction.biomass_components:
            errors.append("Error: submodel '{}' uses dfba but lacks a biomass reaction".format(submodel.name))

        else:
            submodel_species_ids = set([s.id() for s in submodel.get_species()])
            for biomass_component in submodel.biomass_reaction.biomass_components:
                species_id = Species.gen_id(biomass_component.species_type.id,
                                            submodel.biomass_reaction.compartment.id)
                if species_id not in submodel_species_ids:
                    errors.append("Error: undefined species '{}' in biomass reaction '{}' used by "
                                  "submodel '{}'.".format(species_id, submodel.biomass_reaction.name, submodel.name))

        return errors

    def check_dynamic_submodel(self, submodel):
        """ Check the inputs to a dynamic submodel

        Ensure that:

            * All reactions have rate laws for the appropriate directions

        Args:
            submodel (`Submodel`): a dynamic (SSA or ODE) submodel

        Returns:
            :obj:`list` of :obj:`str`: if no errors, returns an empty `list`; otherwise a `list` of
            error messages
        """
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
        """ Transcode and evaluate all rate law equations in a model

        Ensure that all rate law equations can be transcoded and evaluated. Rate laws that
        succesfully transcode are stored in `rate_law.equation.transcoded`.

        Returns:
            :obj:`list` of `str`: if no errors, returns an empty `list`; otherwise a `list` of
            error messages
        """
        errors = []

        species_concentrations = {}
        for concentration in self.model.get_concentrations():
            species_concentrations[concentration.species.serialize()] = concentration.value

        parameters = self.model.get_parameters()
        parameter_ids = set([parameter.id for parameter in parameters])
        parameter_values = {}
        for parameter in parameters:
            parameter_values[parameter.id] = parameter.value

        species_ids = set([specie.id() for specie in self.model.get_species()])
        for reaction in self.model.get_reactions():
            for rate_law in reaction.rate_laws:
                if getattr(rate_law, 'equation', None) is None:
                    continue
                try:
                    rate_law.equation.transcoded = RateLawUtils.transcode(rate_law.equation, species_ids, parameter_ids)
                except Exception as error:
                    errors.append('{} rate law for reaction "{}" cannot be transcoded: {}'.format(
                        rate_law.direction.name, reaction.id, str(error)))
            try:
                rates = RateLawUtils.eval_reaction_rate_laws(reaction, species_concentrations, parameter_values)
            except Exception as error:
                errors.append('{} for reaction "{}" cannot be evaluated: {}'.format(
                    rate_law.direction.name, reaction.id, str(error)))
        return errors

    # TODO(Arthur): reconsider; not good for dFBA models; perhaps good for dynamic models
    def verify_reactant_compartments(self):
        """ Verify that all reactants in each submodel's reactions are in the submodel's compartment

        Returns:
            :obj:`list` of `str`: if no errors, returns an empty `list`; otherwise a `list` of
            error messages
        """
        errors = []
        for submodel in self.model.get_submodels():
            compartment = submodel.compartment
            if compartment is None:
                errors.append("submodel '{}' must contain a compartment attribute".format(submodel.id))
                continue
            for reaction in submodel.reactions:
                for participant in reaction.participants:
                    if participant.coefficient < 0:     # select reactants
                        if participant.species.compartment != compartment:
                            error = "submodel '{}' models compartment {}, but its reaction {} uses "\
                                "specie {} in another compartment: {}".format(
                                    submodel.id,
                                    compartment.id, reaction.id, participant.species.species_type.id,
                                    participant.species.compartment.id)
                            errors.append(error)
        return errors

    def verify_species_types(self):
        """ Verify all species types

        Ensure that:

            * All species types have positive molecular weights

        Returns:
            :obj:`list` of `str`: if no errors, returns an empty `list`; otherwise a `list` of
            error messages
        """
        errors = []
        for species_type in self.model.get_species_types():
            if not 0<species_type.molecular_weight:
                errors.append("species types must contain positive molecular weights, but the MW for {} "
                    "is {}".format(species_type.id, species_type.molecular_weight))
        return []

    def verify_acyclic_dependencies(self, model_types):
        """ Verify that the network of depencencies for model types in `model_types` are acyclic

        Ensure that:
            * The model types in `model_types` do not make recursive calls; tested types
                include Observable and FunctionExpression

        Args:
            model_types (:obj:`list` of `type`): model types (subclasses of `obj_model.Model`) to test

        Returns:
            :obj:`list` of :obj:`str`: if no errors, returns an empty `list`; otherwise a `list` of
            error messages
        """
        errors = []
        for model_type in model_types:

            # get all instances of model_type in self.model
            all_models = None
            for name, attr in self.model.Meta.related_attributes.items():
                if hasattr(self.model.Meta.related_attributes[name], 'primary_class') and \
                    self.model.Meta.related_attributes[name].primary_class == model_type:
                    all_models = getattr(self.model, name)

            # get self-referential attribute, if any
            expression_model = model_type.Meta.expression_model
            name_self_ref_attr = None
            for name in expression_model.Meta.attributes.keys():
                if hasattr(expression_model.Meta.attributes[name], 'related_class') and \
                    expression_model.Meta.attributes[name].related_class == model_type:
                    name_self_ref_attr = name

            if all_models and name_self_ref_attr:
                digraph = nx.DiGraph()
                for model in all_models:
                    for ref in getattr(model.expression, name_self_ref_attr):
                        digraph.add_edge(model, ref)
                cycle_generator = nx.simple_cycles(digraph)
                for cycle in cycle_generator:
                    cyc = [o.id for o in cycle]
                    cyc.append(cyc[0])
                    errors.append("dependency cycle among {}s: {}".format(model_type.__name__, '->'.join(cyc)))
        return errors
