import pkg_resources

with open(pkg_resources.resource_filename('wc_lang', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import (RateLawDirection, TaxonRank,
                   Model, Taxon, Submodel, Compartment, InitVolume, Ph,
                   ChemicalStructure, ChemicalStructureFormat, ChemicalStructureAlphabet,
                   SpeciesType, Species, DistributionInitConcentration, DfbaObjective, DfbaObjectiveExpression,
                   Observable, ObservableExpression,
                   Function, FunctionExpression,
                   Reaction, FluxBounds, SpeciesCoefficient, RateLaw, RateLawExpression,
                   DfbaObjSpecies, DfbaObjReaction, Parameter,
                   StopCondition, StopConditionExpression,
                   Observation, ObservationGenotype, ObservationEnv, ObservationSet, Evidence, Conclusion, Process, 
                   Identifier, Reference, Author, Change,
                   Validator, WcLangWarning)
