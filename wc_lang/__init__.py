import pkg_resources

with open(pkg_resources.resource_filename('wc_lang', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import (TaxonRank, SubmodelAlgorithm, SpeciesTypeType, RateLawDirection, ReferenceType,
                   Model, Taxon, Submodel, ObjectiveFunction, Compartment,
                   SpeciesType, Species, Observable, Concentration,
                   Reaction, SpeciesCoefficient, RateLaw, RateLawEquation,
                   BiomassComponent, BiomassReaction, Parameter, Reference,
                   DatabaseReference)
from . import config
from . import io
from . import model_gen
from . import prepare
from . import rate_law_utils
from . import sbml
from . import transform
from . import util
