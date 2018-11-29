import pkg_resources

with open(pkg_resources.resource_filename('wc_lang', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import (TaxonRank, SubmodelAlgorithm, SpeciesTypeType, ConcentrationUnit, RateLawDirection, ReferenceType,
                   Model, Taxon, Submodel, Compartment,
                   SpeciesType, Species, DfbaObjective, DfbaObjectiveExpression,
                   Observable, ObservableExpression,
                   Function, FunctionExpression, Concentration,
                   Reaction, SpeciesCoefficient, RateLaw, RateLawExpression, Expression,
                   DfbaNetComponent, DfbaNetReaction, Parameter, StopCondition, StopConditionExpression,
                   Reference, DatabaseReference)
from . import config
from . import expression
from . import io
from . import prep_for_sim
from . import sbml
from . import transform
from . import util
