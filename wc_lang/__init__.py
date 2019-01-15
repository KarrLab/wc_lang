import pkg_resources

with open(pkg_resources.resource_filename('wc_lang', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import (Unit, TimeUnit, MassUnit, VolumeUnit, DensityUnit, MoleculeCountUnit, ConcentrationUnit,
                   ReactionParticipantUnit, ReactionRateUnit, ReactionFluxBoundUnit, DfbaObjectiveUnit,
                   DfbaCellSizeUnit, DfbaObjectiveCoefficientUnit, DfbaObjSpeciesUnit,

                   RateLawDirection, TaxonRank,

                   Model, Taxon, Submodel, Compartment,
                   SpeciesType, Species, DistributionInitConcentration, DfbaObjective, DfbaObjectiveExpression,
                   Observable, ObservableExpression,
                   Function, FunctionExpression,
                   Reaction, SpeciesCoefficient, RateLaw, RateLawExpression,
                   DfbaObjSpecies, DfbaObjReaction, Parameter,
                   StopCondition, StopConditionExpression, StopConditionUnit,
                   Evidence, Interpretation, DatabaseReference, Reference,

                   Validator)
from . import config
from . import io
from . import sbml
from . import transform
from . import util
