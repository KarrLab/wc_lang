import pkg_resources

with open(pkg_resources.resource_filename('wc_lang', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import (TaxonRank, SubmodelAlgorithm, SpeciesTypeType, ConcentrationUnit, RateLawDirection, ReferenceType,
                   Model, Taxon, Submodel, ObjectiveFunction, Compartment,
                   SpeciesType, Species, Observable, Function, Concentration,
                   Reaction, SpeciesCoefficient, ObservableCoefficient, RateLaw, RateLawEquation,
                   BiomassComponent, BiomassReaction, Parameter, StopCondition, Reference,
                   DatabaseReference)
from . import config
from . import io
from . import prepare
from . import expression_utils
from . import sbml
from . import transform
from . import util
