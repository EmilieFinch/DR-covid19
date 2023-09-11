# -------------------------------------------------------------------------------------------------------
# Run counterfactual scenarios for Dominican Republic
# Author: Emilie Finch
# -------------------------------------------------------------------------------------------------------

# Run scenarios

no_vacc_counterfactual <- run_counterfactual(base_run = mod_run, sim_end = "2021-12-15", vacc_scenario = "no_vaccination",
                                             REP_START, REP_END, BURN_IN, 
                                             BURN_IN_FINAL, ITER, ITER_FINAL, 
                                             replic, which_pops, params,
                                             ld, sitreps, sero, 
                                             virus, sgtf, set_id, 
                                             opt_mobility, seas_yn, seas_amp, 
                                             opt_v2, opt_relu, voct, voci, vacc)


pfizer_counterfactual <-  run_counterfactual(run_counterfactual(base_run = mod_run, sim_end = "2021-12-15", vacc_scenario = "pfizer-instead",
                                                                REP_START, REP_END, BURN_IN, 
                                                                BURN_IN_FINAL, ITER, ITER_FINAL, 
                                                                replic, which_pops, params,
                                                                ld, sitreps, sero, 
                                                                virus, sgtf, set_id, 
                                                                opt_mobility, seas_yn, seas_amp, 
                                                                opt_v2, opt_relu, voct, voci, vacc))

pfizer_delay_counterfactual <-  run_counterfactual(base_run = mod_run, sim_end = "2021-12-15", vacc_scenario = "pfizer-delay",
                                                   REP_START, REP_END, BURN_IN, 
                                                   BURN_IN_FINAL, ITER, ITER_FINAL, 
                                                   replic, which_pops, params,
                                                   ld, sitreps, sero, 
                                                   virus, sgtf, set_id, 
                                                   opt_mobility, seas_yn, seas_amp, 
                                                   opt_v2, opt_relu, voct, voci, vacc)

az_counterfactual <-  run_counterfactual(base_run = mod_run, sim_end = "2021-12-15", vacc_scenario = "az-instead",
                                         REP_START, REP_END, BURN_IN, 
                                         BURN_IN_FINAL, ITER, ITER_FINAL, 
                                         replic, which_pops, params,
                                         ld, sitreps, sero, 
                                         virus, sgtf, set_id, 
                                         opt_mobility, seas_yn, seas_amp, 
                                         opt_v2, opt_relu, voct, voci, vacc)

az_delay_counterfactual <-  run_counterfactual(base_run = mod_run, sim_end = "2021-12-15", vacc_scenario = "az-delay",
                                               REP_START, REP_END, BURN_IN, 
                                               BURN_IN_FINAL, ITER, ITER_FINAL, 
                                               replic, which_pops, params,
                                               ld, sitreps, sero, 
                                               virus, sgtf, set_id, 
                                               opt_mobility, seas_yn, seas_amp, 
                                               opt_v2, opt_relu, voct, voci, vacc)

scenarios_out <- list(no_vacc = no_vacc_counterfactual,
                           delayed_vacc = delayed_vacc_counterfactual,
                           pfizer = pfizer_counterfactual,
                           pfizer_delay = pfizer_delay_counterfactual,
                           az = az_counterfactual,
                           az_delay = az_delay_counterfactual)
  
qsave(scenarios_out, file = here("output", paste0("vacc-counterfactual_", set_id, ".qs"))) # Save locally


