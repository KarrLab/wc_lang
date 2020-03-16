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
    {'path': 'wc_sim/tests/fixtures/test_model.xlsx',
     'ignore_extra_models': True},
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
    {'path': 'wc_sim/tests/fixtures/verification/cases/multialgorithmic/00001/00001-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/multialgorithmic/00003/00003-wc_lang_1_submodel.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/multialgorithmic/00003/00003-wc_lang_2_submodels.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/multialgorithmic/00003/00003-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/multialgorithmic/00007/00007-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/multialgorithmic/00020/00020-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/multialgorithmic/00021/00021-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/multialgorithmic/00030/00007-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/multialgorithmic/00030/00030-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00001/00001-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00002/00002-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00003/00003-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00004/00004-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00005/00005-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00006/00006-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00010/00010-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00014/00014-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00015/00015-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00017/00017-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00018/00018-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00019/00019-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00020/00020-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00021/00021-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00022/00022-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00028/00028-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/semantic/00054/00054-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00001/00001-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00003/00003-wc_lang_1_submodel.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00003/00003-wc_lang_2_submodels.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00003/00003-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00004/00004-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00007/00007-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00007_hybrid/00007_hybrid-wc_lang_old.xlsx',
     'validate': False},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00007_hybrid/00007_hybrid-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00012/00012-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00020/00020-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00021/00021-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00030/00030-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/00037/00037-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/cases/stochastic/transcription_translation/transcription_translation-wc_lang.xlsx',
     'validate': False},
    {'path': 'wc_sim/tests/fixtures/verification/testing/hybrid/transcription_translation/transcription_translation_correct_ssa.xlsx',
     'validate': False},
    {'path': 'wc_sim/tests/fixtures/verification/testing/hybrid/transcription_translation/transcription_translation_hybrid.xlsx',
     'validate': False},
    {'path': 'wc_sim/tests/fixtures/verification/testing/hybrid/transcription_translation/transcription_translation-wc_lang_JK.xlsx',
     'validate': False},
    {'path': 'wc_sim/tests/fixtures/verification/testing/hybrid/translation_metabolism/translation_metabolism_correct_ssa.xlsx',
     'validate': False},
    {'path': 'wc_sim/tests/fixtures/verification/testing/hybrid/translation_metabolism/translation_metabolism_hybrid.xlsx',
     'validate': False},
    {'path': 'wc_sim/tests/fixtures/verification/testing/multialgorithmic/00007/00007-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/testing/semantic/00001/00001-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/testing/semantic/00004/00004-wc_lang.xlsx',
     'validate': False},
    {'path': 'wc_sim/tests/fixtures/verification/testing/semantic/00054/00054-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/testing/stochastic/00001/00001-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/testing/stochastic/00006/00006-wc_lang.xlsx'},
    {'path': 'wc_sim/tests/fixtures/verification/testing_ValidationSuite_run/stochastic/00001/00001-wc_lang.xlsx',
     'validate': False},
    {'path': 'wc_sim/tests/fixtures/verification/testing_ValidationSuite_run/stochastic/00006/00006-wc_lang.xlsx',
     'validate': False},

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
    {'path': 'h1_hesc/tests/code/fixtures/mock_model.xlsx'},

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
    if path.get('validate', True):
        kwargs = copy.copy(path)
        kwargs.pop('path')
        if 'validate' in kwargs:
            kwargs.pop('validate')
        try:
            wc_lang.io.Reader().run(abs_path, **kwargs)
        except ValueError as err:
            warnings.warn('{} is invalid: {}'.format(path['path'], str(err)))
