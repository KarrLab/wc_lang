import migration_2020_03_09 as migration
""" Migration WC-lang-encoded files

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-10-10
:Copyright: 2019, Karr Lab
:License: MIT
"""

import copy
import os.path
import sys
import warnings
import wc_lang.io
sys.path.insert(0, 'migrations')

base_dir = os.path.expanduser('~/Documents')

paths = [
    # wc_lang
    {'path': 'wc_lang/tests/fixtures/example-model.xlsx'},
    {'path': 'wc_lang/tests/fixtures/sbml-io.xlsx'},
    {'path': 'wc_lang/tests/fixtures/sbml-io-transformed.xlsx'},
    {'path': 'wc_lang/tests/fixtures/test_model.xlsx'},
    {'path': 'wc_lang/tests/fixtures/test_validate_model.xlsx'},
    {'path': 'wc_lang/tests/sbml/fixtures/static-model.xlsx'},

    # wc_sim
    {'path': 'wc_sim/examples/transcription_translation_hybrid_model/model.xlsx'},
    {'path': 'wc_sim/examples/translation_metabolism_hybrid_model/model.xlsx'},
    {'path': 'wc_sim/tests/fixtures/2_species_1_reaction.xlsx'},
    {'path': 'wc_sim/tests/fixtures/2_species_1_reaction_with_rates_given_by_reactant_population.xlsx'},
    {'path': 'wc_sim/tests/fixtures/2_species_a_pair_of_symmetrical_reactions_rates_given_by_reactant_population.xlsx'},
    {'path': 'wc_sim/tests/fixtures/MetabolismAndGeneExpression.xlsx'},
    {'path': 'wc_sim/tests/fixtures/test_dry_model.xlsx'},
    {'path': 'wc_sim/tests/fixtures/test_dry_model_with_mass_computation.xlsx',
     'ignore_extra_models': True},
    {'path': 'wc_sim/tests/fixtures/test_dynamic_expressions.xlsx'},
    {'path': 'wc_sim/tests/fixtures/test_model.xlsx'},
    {'path': 'wc_sim/tests/fixtures/test_model_for_access_species_populations.xlsx'},
    {'path': 'wc_sim/tests/fixtures/test_model_for_access_species_populations_steady_state.xlsx'},
    {'path': 'wc_sim/tests/fixtures/test_new_features_model.xlsx'},
    {'path': 'wc_sim/tests/fixtures/dynamic_tests/one_exchange_rxn_compt_growth.xlsx',
     'ignore_extra_models': True},
    {'path': 'wc_sim/tests/fixtures/dynamic_tests/stop_conditions.xlsx',
     'ignore_extra_models': True},
    {'path': 'wc_sim/tests/fixtures/dynamic_tests/one_reaction_linear.xlsx',
     'ignore_extra_models': True},
    {'path': 'wc_sim/tests/fixtures/dynamic_tests/template.xlsx',
     'ignore_extra_models': True},
    {'path': 'wc_sim/tests/fixtures/dynamic_tests/one_rxn_exponential.xlsx',
     'ignore_extra_models': True},
    {'path': 'wc_sim/tests/fixtures/dynamic_tests/static.xlsx',
     'ignore_extra_models': True},
    {'path': 'wc_sim/tests/submodels/fixtures/test_submodel.xlsx'},
    {'path': 'wc_sim/tests/submodels/fixtures/test_submodel_no_shared_species.xlsx'},

    # wc_test
    {'path': 'wc_test/tests/fixtures/min_model.xlsx'},

    # intro_to_wc_modeling
    {'path': 'intro_to_wc_modeling/intro_to_wc_modeling/wc_modeling/wc_lang_tutorial/examples/example_model.xlsx'},

    # wc_analysis
    {'path': 'wc_analysis/tests/fixtures/test_model.xlsx'},

    # mycoplasma_pneumoniae
    {'path': 'mycoplasma_pneumoniae/mycoplasma_pneumoniae/model/model_calibration.xlsx'},
    {'path': 'mycoplasma_pneumoniae/mycoplasma_pneumoniae/model/model_calibration_wDeg.xlsx'},

    # h1_hesc
    {'path': 'h1_hesc/tests/model/cell_cycle/fixtures/test_exponential_growth_in_M.xlsx'},
    {'path': 'h1_hesc/tests/model/cell_cycle/fixtures/test_cyclin_dynamics.xlsx'},
    {'path': 'h1_hesc/model/hesc_recon/hesc_recon_wc_data_model.xlsx'},
    {'path': 'h1_hesc/tests/code/fixtures/eukaryote_model.xlsx'},

    # rand_wc_model_gen
    {'path': 'rand_wc_model_gen/rand_wc_model_gen/model_gen/model.xlsx'},
    {'path': 'rand_wc_model_gen/rand_wc_model_gen/model_gen/model_2.xlsx'},
]

for i_path, path in enumerate(paths):
    print('Migrating path {} of {}: {}'.format(i_path + 1, len(paths), path['path']))

    abs_path = os.path.join(base_dir, path['path'])

    # migrate
    migration.transform(abs_path)

    # validate
    kwargs = copy.copy(path)
    kwargs.pop('path')
    try:
        wc_lang.io.Reader().run(abs_path, **kwargs)
    except ValueError as err:
        warnings.warn('{} is invalid: {}'.format(path['path'], str(err)))
