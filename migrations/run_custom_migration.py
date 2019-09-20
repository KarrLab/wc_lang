cmd='et'
cmd='python3.6 migrations/2019_09_20.py'
cmd='wc-lang validate'

# wc_lang
    ${cmd} ~/Documents/wc_lang/tests/fixtures/example-model.xlsx
    ${cmd} ~/Documents/wc_lang/tests/fixtures/test_model.xlsx
    ${cmd} ~/Documents/wc_lang/tests/fixtures/test_validate_model.xlsx
    ${cmd} ~/Documents/wc_lang/tests/fixtures/sbml-io.xlsx
    ${cmd} ~/Documents/wc_lang/tests/fixtures/sbml-io-transformed.xlsx

# wc_sim
    ${cmd} ~/Documents/wc_sim/examples/transcription_translation_hybrid_model/model.xlsx
    ${cmd} ~/Documents/wc_sim/examples/transcription_translation_hybrid_model/model_for_migration.bak.xlsx
    ${cmd} ~/Documents/wc_sim/examples/transcription_translation_hybrid_model/model_for_migration.xlsx
    # ${cmd} ~/Documents/wc_sim/examples/transcription_translation_hybrid_model/model_migrated_forward.xlsx
    # ${cmd} ~/Documents/wc_sim/examples/transcription_translation_hybrid_model/model_migrated_roundtrip.xlsx
    ${cmd} ~/Documents/wc_sim/examples/translation_metabolism_hybrid_model/model.xlsx
    ${cmd} ~/Documents/wc_sim/tests/fixtures/2_species_1_reaction.xlsx
    ${cmd} ~/Documents/wc_sim/tests/fixtures/2_species_1_reaction_with_rates_given_by_reactant_population.xlsx
    ${cmd} ~/Documents/wc_sim/tests/fixtures/2_species_a_pair_of_symmetrical_reactions_rates_given_by_reactant_population.xlsx
    ${cmd} ~/Documents/wc_sim/tests/fixtures/MetabolismAndGeneExpression.xlsx
    ${cmd} ~/Documents/wc_sim/tests/fixtures/test_dry_model.xlsx
    # ${cmd} ~/Documents/wc_sim/tests/fixtures/test_dry_model_with_mass_computation.xlsx
    ${cmd} ~/Documents/wc_sim/tests/fixtures/test_model.xlsx
    ${cmd} ~/Documents/wc_sim/tests/fixtures/test_model_for_access_species_populations.xlsx
    ${cmd} ~/Documents/wc_sim/tests/fixtures/test_model_for_access_species_populations_steady_state.xlsx
    # ${cmd} ~/Documents/wc_sim/tests/fixtures/test_model_with_mass_computation.xlsx
    ${cmd} ~/Documents/wc_sim/tests/fixtures/test_new_features_model.xlsx
    ${cmd} ~/Documents/wc_sim/tests/submodels/fixtures/test_submodel.xlsx
    ${cmd} ~/Documents/wc_sim/tests/submodels/fixtures/test_submodel_no_shared_species.xlsx

# wc_test
    ${cmd} ~/Documents/wc_test/tests/fixtures/min_model.xlsx

# intro_to_wc_modeling
    ${cmd} ~/Documents/intro_to_wc_modeling/intro_to_wc_modeling/wc_modeling/wc_lang_tutorial/examples/example_model.xlsx

# wc_analysis
    ${cmd} ~/Documents/wc_analysis/tests/fixtures/test_model.xlsx

# mycoplasma_pneumoniae
    ${cmd} ~/Documents/mycoplasma_pneumoniae/mycoplasma_pneumoniae/model/model_calibration.xlsx
    ${cmd} ~/Documents/mycoplasma_pneumoniae/mycoplasma_pneumoniae/model/model_calibration_wDeg.xlsx

# h1_hesc
    ${cmd} ~/Documents/h1_hesc/tests/model/cell_cycle/fixtures/test_exponential_growth_in_M.xlsx
    ${cmd} ~/Documents/h1_hesc/tests/model/cell_cycle/fixtures/test_cyclin_dynamics.xlsx
    CONFIG__DOT__wc_lang__DOT__validation__DOT__validate_element_charge_balance=0 ${cmd} ~/Documents/h1_hesc/model/hesc_recon/hesc_recon_wc_data_model.xlsx

# rand_wc_model_gen
    ${cmd} ~/Documents/rand_wc_model_gen/rand_wc_model_gen/model_gen/model.xlsx
    ${cmd} ~/Documents/rand_wc_model_gen/rand_wc_model_gen/model_gen/model_2.xlsx
    