# Function to fit covidM to data

## This function takes fits the compartmental model (specified in sim_compartment) to data using MCMC.
## Input parameters include parameters defining burn-in and iterations of the MCMC chain, priors, model parameters and data to fit to 
## The function returns test (a summary of dynamicscover time by age, time and run), posteriorsI (containing posterior values for each trial), parametersI, priorsI, and constants (used in the model fitting).

covidm_fit <- function(REP_START, REP_END, BURN_IN, 
                       BURN_IN_FINAL, ITER, ITER_FINAL, 
                       replic, which_pops, priorsI,
                       params, ld, sitreps, sero, 
                       virus, sgtf, set_id, opt_mobility, 
                       seas_yn, seas_amp, opt_v2, 
                       opt_relu, voct, voci, vacc, date_fitting){
  
  posteriorsI <- list()
  dynamicsI <- list()
  parametersI <- list()
  

  # Check execution time
  start_time <- Sys.time()
  
  # Loop through regions
  for (pn in which_pops) {
    paramsI <- rlang::duplicate(params)
    paramsI$pop <- list(rlang::duplicate(params$pop[[pn]]))
    paramsI$travel <- matrix(1, nrow = 1, ncol = 1)
    paramsI$schedule <- list()
    j = 1
    
    for (i in seq_along(params$schedule)) {
      if (pn - 1 == params$schedule[[i]]$pops) {
        paramsI$schedule[[j]] <- rlang::duplicate(params$schedule[[i]])
        paramsI$schedule[[j]]$pops <- 0
        j = j + 1
      }
    }
    
    ldI <- rlang::duplicate(ld)
    ldI <- ldI[pid == pn - 1]
    sitrepsI <- rlang::duplicate(sitreps)
    sitrepsI <- sitrepsI[pid == pn - 1]
    seroI <- rlang::duplicate(sero)
    seroI <- seroI[pid == pn - 1]   
    virusI <- rlang::duplicate(virus)
    virusI <- virusI[pid == pn - 1]
    sgtfI <- copy(sgtf)
    sgtfI <- sgtfI[pid == pn - 1]
    
    # load cpp functions
    cm_source_backend(
      user_defined = list(
        model_v2 = list(
          
          cpp_changes = cpp_chgI_voc(priorsI, constants, 
                                     v2 = opt_v2, v2_relu = opt_relu, mobility_type = opt_mobility),
          
          cpp_loglikelihood = cpp_likI_voc(paramsI, ldI, sitrepsI, 
                                           virusI, sgtfI, seroI, 
                                           pn, date_fitting, priorsI, 
                                           constants, death_cutoff = 0, use_sgtf = FALSE),
          
          cpp_observer = c(cpp_obsI_voc(v2 = opt_v2, P.death, P.critical, 
                                        P.hosp, priorsI, constants, voct, voci),
                           cpp_obsI_vax(paramsI, vacc[[pn]]),
                           if(seas_yn == "seasyes") cpp_obsI_seasonality(seas_amp,1) else "")
          
        )))
    
    priorsI2 <- rlang::duplicate(priorsI)

    postI <- cm_backend_mcmc_test(cm_translate_parameters(paramsI), priorsI2,
                                 seed = 0, 
                                 burn_in = ifelse(replic == REP_END, BURN_IN_FINAL, BURN_IN), 
                                 iterations = ifelse(replic == REP_END, ITER_FINAL, ITER), 
                                 n_threads = N_THREADS, classic_gamma = F, do_migration = T)
    setDT(postI)
    postI = cbind(postI, as.data.table(constants))
    posteriorsI[[pn]] = postI
    parametersI[[pn]] = rlang::duplicate(paramsI)
    print(pn)
  }
  
  # Check time
  time2 <- Sys.time()
  print(time2-start_time)

  # Sample dynamics from fit
  dynamicsI <- list()
  for (pn in which_pops)  {
    cat(paste0("Sampling fit for population ", pn, "...\n"))
    
    # Source backend
    cm_source_backend(
      user_defined = list(
        model_v2 = list(
          
          cpp_changes = cpp_chgI_voc(priorsI, constants, 
                                     v2 = opt_v2, v2_relu = opt_relu, mobility_type = opt_mobility),
          
          cpp_loglikelihood = "",
          
          cpp_observer = c(cpp_obsI_voc(v2 = opt_v2, P.death, P.critical, 
                                        P.hosp, priorsI, constants, voct, voci),
                           cpp_obsI_vax(parametersI[[pn]], vacc[[pn]]),
                           if (seas_yn == "seasyes") cpp_obsI_seasonality(seas_amp, 1) else "")
        )))
    
    # Sampling fits
    paramsI2 <- rlang::duplicate(parametersI[[pn]])
    paramsI2$time1 <- as.character(ymd(parametersI[[pn]]$time1) + 56)
    test <- cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriorsI[[pn]], 500, seed = 0)
    rows <- cm_backend_sample_fit_rows(cm_translate_parameters(paramsI2), posteriorsI[[pn]], 500, seed = 0)
    
    test <- rbindlist(test)
    test[, population := pn]
    
    # Add dispersion parameters
    disp = posteriorsI[[pn]][rows, .SD, .SDcols = patterns("^disp|v2_conc|v2_disp|v2_sgtf0")]
    disp[, run := .I]
    test = merge(test, disp, by = "run")
    
    dynamicsI[[pn]] = test
    
  }
  
  # Concatenate dynamics 
  
  test = rbindlist(dynamicsI, fill = TRUE)
  test$population <- as.character(test$population)
  test[, population := "Dominican Republic"]
  qsave(rlang::duplicate(list(posteriorsI, parametersI, priorsI, test, constants)), here("output",paste0(set_id, ".qs")))

  return(list(dynamics = test, posteriors = posteriorsI, parameters = parametersI, priors = priorsI, constants = constants))
 
}
  