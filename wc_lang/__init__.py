import pkg_resources

with open(pkg_resources.resource_filename('wc_lang', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from .core import (TimeUnit, TaxonRank,
                   CompartmentBiologicalType, CompartmentPhysicalType, CompartmentGeometry,
                   MassUnit, VolumeUnit, DensityUnit,
                   SubmodelAlgorithm, SpeciesTypeType, MoleculeCountUnit,
                   RandomDistribution, ConcentrationUnit,
                   ReactionParticipantUnit, RateLawDirection, RateLawType, ReactionRateUnit, ReactionFluxBoundUnit,
                   DfbaObjectiveUnit, DfbaCellSizeUnit, DfbaObjectiveCoefficientUnit, DfbaNetComponentUnit,
                   ParameterType, EvidenceType, ReferenceType,
                   Model, Taxon, Submodel, Compartment,
                   SpeciesType, Species, DistributionInitConcentration, DfbaObjective, DfbaObjectiveExpression,
                   Observable, ObservableExpression,
                   Function, FunctionExpression,
                   Reaction, SpeciesCoefficient, RateLaw, RateLawExpression,
                   DfbaNetSpecies, DfbaNetReaction, Parameter,
                   StopCondition, StopConditionExpression, StopConditionUnit,
                   Evidence, DatabaseReference, Reference,
                   Validator)
from . import config
from . import io
from . import sbml
from . import transform
from . import util
