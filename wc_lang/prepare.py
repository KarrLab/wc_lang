""" Prepare a WC model for further processing, such as export or simulation

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-10-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

from math import ceil, floor, exp, log, log10, isnan
from obj_model import utils
from wc_utils.util.list import difference
from wc_lang.core import (SubmodelAlgorithm, Model, ObjectiveFunction, SpeciesType, SpeciesTypeType,
    Species, Compartment, Reaction, ReactionParticipant, RateLawEquation, BiomassReaction)
from wc_sim.multialgorithm.model_utilities import ModelUtilities

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
        Fill gaps in dFBA submodel reaction networks
        Ensure that dFBA submodels have objective functions
        Apply default flux bounds to the reactions in dFBA submodels
    '''
    def __init__(self, model):
        self.model = model

    def run(self):
        for submodel in self.model.get_submodels():
            if submodel.algorithm == SubmodelAlgorithm.dfba:
                reactions_created = self.fill_dfba_submodel_reaction_gaps(submodel)
                self.confirm_dfba_submodel_obj_func(submodel)
                (min_bounds_set, min_bounds_set) = self.apply_default_dfba_submodel_flux_bounds(submodel)

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
            ValueError: if some species in `submodel` are neither produced nor consumed

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
            raise ValueError("some species in submodel '{}' are neither produced nor consumed: {}".format(
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
            return

        # use the biomass reaction
        obj_func_expression = submodel.biomass_reaction.id
        # deserialize the expression
        attr = ObjectiveFunction.Meta.attributes['expression']
        # deserialize needs the biomass reaction and all the Reactions
        objs = {}
        objs[BiomassReaction] = {submodel.biomass_reaction.id:submodel.biomass_reaction}
        objs[Reaction] = dict(zip([rxn.id for rxn in submodel.reactions], submodel.reactions))
        (of, invalid_attribute) = ObjectiveFunction.deserialize(attr, obj_func_expression, objs)
        if invalid_attribute:
            raise ValueError("submodel '{}' cannot use biomass reaction '{}' as an objective function: "
                "{}".format(submodel.name, submodel.biomass_reaction.id, invalid_attribute.messages[0]))
        submodel.objective_function = of

    # TODO: doublecheck Sphinx formatting
    def apply_default_dfba_submodel_flux_bounds(self, submodel):
        ''' Apply default flux bounds to a dFBA submodel's reactions

        The FBA optimizer needs min and max flux bounds for each dFBA submodel reaction.
        If bounds are not provided in some reactions, and default bounds are provided in a config file,
        then apply the defaults to the reactions.
        Specifically, min and max default bounds are applied as follows:
            reversible reactions:
                min_flux = -default_max_flux_bound
                max_flux = default_max_flux_bound
            irreversible reactions:
                min_flux = default_min_flux_bound
                max_flux = default_max_flux_bound

        Args:
            submodel (`Submodel`): a dFBA submodel

        Raises:
            ValueError: if `submodel` is not a dFBA submodel

        Returns:
            :obj:`tuple` of (`int`,`int`): counts of min and max flux bounds set
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
            * The DFBA submodel contains a biomass reaction

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
        # TODO: avoid this redundancy:
        This method is deliberately redundant with `MultialgorithmSimulation.transcode_rate_laws()`,
        which does not report errors.

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
                    rate_law.equation.transcoded = ModelUtilities.transcode(rate_law, species)
                except Exception as error:
                    errors.append(str(error))
            try:
                rates = ModelUtilities.eval_rate_laws(reaction, species_concentrations)
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
