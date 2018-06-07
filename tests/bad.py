import os

from wc_lang.io import Reader, Writer
from wc_lang.core import (Model, Submodel,  SpeciesType, SpeciesTypeType, Species,
                          Reaction, Observable, Compartment,
                          SpeciesCoefficient, ObservableCoefficient, Parameter,
                          RateLaw, RateLawDirection, RateLawEquation, SubmodelAlgorithm, Concentration,
                          BiomassComponent, BiomassReaction, StopCondition, ConcentrationUnit)

model = Model(id='test_model', version='0.0.0')
species_type = model.species_types.create(id='st')
comp = model.compartments.create(id='compt')
spec = comp.species.create(species_type=species_type)
# the round trip model is the same if coeff != 1
coeff = 1
species_coefficient = spec.species_coefficients.create(coefficient=coeff)
obs_plain = model.observables.create(id='obs_1')
obs_plain.species.append(species_coefficient) # IMHO should be obs_plain.species_coefficients.append ...
submodel = model.submodels.create(id='s', compartment=comp)
reaction = submodel.reactions.create(id='r')
reaction.participants.create(species=spec, coefficient=1)

assert model.validate() is None

filename = os.path.join('tmp', 'test_model.xlsx')
Writer().run(model, filename, set_repo_metadata_from_path=False)
round_trip_model = Reader().run(filename)
assert round_trip_model.validate() is None
# assert round_trip_model.is_equal(model)
# assert model.difference(round_trip_model) == ''
print('', model.difference(round_trip_model))
