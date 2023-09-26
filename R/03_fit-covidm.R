# -------------------------------------------------------------------------------------------------------
# Fit model for Dominican Republic
# Author: Emilie Finch, adapted from code by Nick Davies in https://github.com/nicholasdavies/newcovid
# -------------------------------------------------------------------------------------------------------

today <- Sys.Date()
set_id <- paste0("fit-", opt_mobility, "_takeoff-", extra_voc_takeoff, "_waning-", opt_wane, "_mob-", opt_mobility, "_", today);

# Set priors for fitted parameters ----------------------------------------------------------------------------
# See Table 3 in manuscript for description and references

priorsI = list(
    u = "N 0.09 0.02 T 0.055 0.2", 
    death_mean = "N 15 2 T 5 30",    
    hosp_admission = "N 8 1 T 4 20", 
    icu_admission = "N 12.5 1 T 8 14", 
    cfr_rlo = "N 0 0.1 T -2 2",
    cfr_rlo2 = "N 0 0.1 T -2 2",
    cfr_rlo3 = "N 0 0.1 T -2 2",
    hosp_rlo = "N 0 0.1 T -2 2", 
    hosp_rlo2 = "N 0 0.1 T -2 2", 
    icu_rlo = "N 0 0.1 T -2 2",
    icu_rlo2 = "N 0 0.1 T -2 2",
    contact_adj_a = "B 15 1", 
    contact_adj_b = "B 15 1", 
    disp_deaths = "E 10 10",
    disp_hosp_inc = "E 10 10",
    disp_hosp_prev = "E 10 10",
    disp_icu_prev = "E 10 10",
    extra_voc_relu = "L 0.4 0.1 T 1 4"
)


if (opt_v2) {
    priorsI = c(priorsI, list(
        v2_relu = "L 1.03 0.01 T 1 4", 
        v2_when = "U  486 517", # mid-April to mid-July  
        v2_disp = "E 10 10 T 0 0.7",
        v2_hosp_rlo = "N 0 0.1 T -4 4",
        v2_icu_rlo = "N 0 0.1 T -4 4",
        v2_cfr_rlo = "N 0 0.1 T -4 4"
    ))
}

constants = list(
  tS = 20)


for (replic in REP_START:REP_END)
{
  mod_run <- covidm_fit(REP_START, REP_END, BURN_IN, 
                    BURN_IN_FINAL, ITER, ITER_FINAL, 
                    replic, which_pops, priorsI, 
                    params, ld, sitreps, sero, 
                    virus, sgtf, set_id, opt_mobility,
                    seas_yn, seas_amp, opt_v2, 
                    opt_relu, voct, voci, vacc, date_fitting)
  
    # Run convergence checks 
    conv_out <- convergence(mod_run$posteriorsI, 1)
    qsave(conv_out, file = here("output", paste0("./conv-diag-", set_id, ".qs"))) # Save locally

    # Generate model fit 
    mod_fit <- gen_fit(mod_run$test, mod_run$parametersI, ld, sitreps, virus, sero, "Dominican Republic")
    qsave((list(posteriors = mod_run$posteriorsI, parameters = mod_run$parametersI, priors = mod_run$priorsI, dynamics = mod_run$test, constants = mod_run$constants, model_fit = mod_fit)), here("output",paste0(set_id, ".qs")))
    
}
