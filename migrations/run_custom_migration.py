""" Migration WC-lang-encoded files

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-10-10
:Copyright: 2019, Karr Lab
:License: MIT
"""

import os.path
import sys
import warnings
import wc_lang.io
sys.path.insert(0, 'migrations')
from migration_2019_10_11 import transform

base_dir = os.path.expanduser('~/Documents')

paths = [
    # wc_lang
    'wc_lang/tests/fixtures/example-model.xlsx',
    'wc_lang/tests/fixtures/test_model.xlsx',
    'wc_lang/tests/fixtures/test_validate_model.xlsx',
    'wc_lang/tests/fixtures/sbml-io.xlsx',
    'wc_lang/tests/fixtures/sbml-io-transformed.xlsx',

    # wc_sim
    'wc_sim/examples/transcription_translation_hybrid_model/model.xlsx',
    'wc_sim/examples/transcription_translation_hybrid_model/model_for_migration.bak.xlsx',
    'wc_sim/examples/transcription_translation_hybrid_model/model_for_migration.xlsx',
    # 'wc_sim/examples/transcription_translation_hybrid_model/model_migrated_forward.xlsx',
    # 'wc_sim/examples/transcription_translation_hybrid_model/model_migrated_roundtrip.xlsx',
    'wc_sim/examples/translation_metabolism_hybrid_model/model.xlsx',
    'wc_sim/tests/fixtures/2_species_1_reaction.xlsx',
    'wc_sim/tests/fixtures/2_species_1_reaction_with_rates_given_by_reactant_population.xlsx',
    'wc_sim/tests/fixtures/2_species_a_pair_of_symmetrical_reactions_rates_given_by_reactant_population.xlsx',
    'wc_sim/tests/fixtures/MetabolismAndGeneExpression.xlsx',
    'wc_sim/tests/fixtures/test_dry_model.xlsx',
    # 'wc_sim/tests/fixtures/test_dry_model_with_mass_computation.xlsx',
    'wc_sim/tests/fixtures/test_model.xlsx',
    'wc_sim/tests/fixtures/test_model_for_access_species_populations.xlsx',
    'wc_sim/tests/fixtures/test_model_for_access_species_populations_steady_state.xlsx',
    # 'wc_sim/tests/fixtures/test_model_with_mass_computation.xlsx',
    'wc_sim/tests/fixtures/test_new_features_model.xlsx',
    'wc_sim/tests/submodels/fixtures/test_submodel.xlsx',
    'wc_sim/tests/submodels/fixtures/test_submodel_no_shared_species.xlsx',

    # wc_test
    'wc_test/tests/fixtures/min_model.xlsx',

    # intro_to_wc_modeling
    'intro_to_wc_modeling/intro_to_wc_modeling/wc_modeling/wc_lang_tutorial/examples/example_model.xlsx',

    # wc_analysis
    'wc_analysis/tests/fixtures/test_model.xlsx',

    # mycoplasma_pneumoniae
    'mycoplasma_pneumoniae/mycoplasma_pneumoniae/model/model_calibration.xlsx',
    'mycoplasma_pneumoniae/mycoplasma_pneumoniae/model/model_calibration_wDeg.xlsx',

    # h1_hesc
    'h1_hesc/tests/model/cell_cycle/fixtures/test_exponential_growth_in_M.xlsx',
    'h1_hesc/tests/model/cell_cycle/fixtures/test_cyclin_dynamics.xlsx',
    'h1_hesc/model/hesc_recon/hesc_recon_wc_data_model.xlsx',
    'h1_hesc/tests/code/fixtures/eukaryote_model.xlsx',

    # rand_wc_model_gen
    'rand_wc_model_gen/rand_wc_model_gen/model_gen/model.xlsx',
    'rand_wc_model_gen/rand_wc_model_gen/model_gen/model_2.xlsx',
]

for i_path, path in enumerate(paths):
    print('Migrating path {} of {}: {}'.format(i_path + 1, len(paths), path))

    abs_path = os.path.join(base_dir, path)

    # migrate
    transform(abs_path)

    # validate
    try:
        wc_lang.io.Reader().run(abs_path)
    except ValueError as err:
        warnings.warn('{} is invalid: {}'.format(path, str(err)))

