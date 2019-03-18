import pkg_resources

with open(pkg_resources.resource_filename('wc_lang', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import (RateLawDirection, TaxonRank,
                   Model, Taxon, Submodel, Compartment,
                   SpeciesType, Species, DistributionInitConcentration, DfbaObjective, DfbaObjectiveExpression,
                   Observable, ObservableExpression,
                   Function, FunctionExpression,
                   Reaction, SpeciesCoefficient, RateLaw, RateLawExpression,
                   DfbaObjSpecies, DfbaObjReaction, Parameter,
                   StopCondition, StopConditionExpression,
                   Evidence, Interpretation, DatabaseReference, Reference, Author, Change,
                   Validator)
