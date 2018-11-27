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
                   BiomassComponent, BiomassReaction, Parameter, StopCondition, StopConditionExpression,
                   Reference, DatabaseReference)
from . import config
from . import io
from . import prepare
from . import expression
from . import sbml
from . import transform
from . import util
