# Functions for counterfactual analysis

# Function to calculate the averted burden for specified variable (e.g. deaths, ICU admissions or hospitalisations)
calc_averted_burden <- function(summ, varname, varname_scen, outname){
  
  quants <- seq(0.05, 0.95, by = 0.05)
  summ[, Baseline := get(varname)]
  summ[, Scenario := get(varname_scen)]
  summ[, Averted := (Scenario - Baseline)]
  qsumm <- melt(summ, id.vars = c("run", "t", "population", "age_group"), measure.vars = c("Baseline", "Scenario", "Averted"))
  tsumm <- qsumm[, as.list(quantile(value, quants, na.rm = TRUE)), by = .(t, population, age_group, variable)]
  
  tsumm[, date := t + ymd("2020-01-01")]
  tsumm <- rbind(tsumm,
                tsumm[, lapply(.SD, sum), .SDcols = `5%`:`95%`, by = .(t, population, date, variable, age_group = rep("All", nrow(tsumm)))],
                use.names = TRUE, fill = TRUE)
  
  data_out <- tsumm[,
                   .(date, `Day of Value` = day(date), `Month of Value` = month(date), `Year of Value` = year(date),
                     AgeBand = age_group, Geography = population, ValueType = outname, ModelType = variable, Value = `50%`,
                     `5%`, `10%`, `15%`, `20%`, `25%`, `30%`, `35%`, `40%`, `45%`, `50%`, `55%`, `60%`, `65%`, `70%`, `75%`, `80%`, `85%`, `90%`, `95%`)]
  
  names(data_out)[names(data_out) %like% "\\%$"] = paste("Quantile", quants);
  
  # names(data_out)[names(data_out) %like% "\\%$"] = paste("Quantile", quants);
  data_out
}

calc_all_averted_burden <- function(btest, stest){
  
  test_all <- merge(btest, stest, by = c("run", "t", "group"), suffixes = c("", "_scen")) 
  
  cat("Adding age groups...\n");
  test_all[group %between% c(1, 1), age_group := "0-4"]
  test_all[group %between% c(2, 3), age_group := "5-14"]
  test_all[group %between% c(4, 5), age_group := "15-24"]
  test_all[group %between% c(6, 9), age_group := "25-44"]
  test_all[group %between% c(10, 13), age_group := "45-64"]
  test_all[group %between% c(14, 15), age_group := "65-74"]
  test_all[group %between% c(16, 16), age_group := "75+"]
  test_all[, age_group := factor(age_group, levels = unique(age_group))]
  
  
  cat("Summarizing variables...\n");
  summ <- test_all[, .(death_o = sum(death_o + death2_o), death_o_scen = sum(death_o_scen + death2_o_scen), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
  summ3 <- test_all[, .(icu_p = sum(icu_p + icu2_p), icu_p_scen = sum(icu_p_scen + icu2_p_scen), disp = mean(disp_icu_prev)), by = .(run, t, population, age_group)]
  summ4 <- test_all[, .(bed_p = sum(pmax(0, hosp_p + hosp2_p - hosp_undetected_p - hosp_undetected2_p)), bed_p_scen = sum(pmax(0, hosp_p_scen + hosp2_p_scen - hosp_undetected_p_scen - hosp_undetected2_p_scen)), disp = mean(disp_hosp_prev)), by = .(run, t, population, age_group)]
  summ34 <- test_all[, .(admissions = sum(hosp_undetected_o + hosp_undetected2_o), admissions_scen = sum(hosp_undetected_o_scen + hosp_undetected2_o_scen), disp = disp_hosp_inc), by = .(run, t, population, age_group)]
  summ_inf_i <- test_all[, .(infections_i = sum(pcr_positive_i), infections_i_scen = sum(pcr_positive_i_scen)), by = .(run, t, population, age_group)]
  summ_inf_p <- test_all[, .(infections_p = sum(pcr_positive_p), infections_p_scen = sum(pcr_positive_p_scen)), by = .(run, t, population, age_group)];
  summ_sero_p <- test_all[, .(sero_p = sum(lfia_positive_p), sero_p_scen = sum(lfia_positive_p_scen)), by = .(run, t, population, age_group)];
  summ_icu_ad = test_all[, .(icu_i = sum(icu_i + icu2_i), icu_i_scen = sum(icu_i_scen + icu2_i_scen)),  by = .(run, t, population, age_group)] 
  summ_hosp_ad = test_all[, .(admissions = sum(hosp_undetected_o + hosp_undetected2_o), admissions_scen = sum(hosp_undetected_o_scen + hosp_undetected2_o_scen), disp = disp_hosp_inc), by = .(run, t, population, age_group)]
  
  cat("Running quantiles...\n");
  w = rbind(
    calc_averted_burden(summ, "death_o", "death_o_scen", "deaths"),
    calc_averted_burden(summ3, "icu_p", "icu_p_scen", "icu_beds"),
    calc_averted_burden(summ4, "bed_p", "bed_p_scen", "hospital_beds"),
    calc_averted_burden(summ34, "admissions", "admissions_scen", "hospital_ad"),
    calc_averted_burden(summ_inf_p, "infections_p", "infections_p", "prevalence_mtp_averted"),
    calc_averted_burden(summ_inf_i, "infections_i", "infections_i", "infections_inc_averted"),
    calc_averted_burden(summ_sero_p, "sero_p", "sero_p_scen", "sero_prev"),
    calc_averted_burden(summ_icu_ad, "icu_i", "icu_i_scen", "icu_i"),
    calc_averted_burden(summ_hosp_ad, "admissions", "admissions_scen", "hosp_ad"))
  
  return (w)
  
}
run_counterfactual <- function(base_run, sim_end, vacc_scenario, 
                              REP_START, REP_END, BURN_IN, 
                              BURN_IN_FINAL, ITER, ITER_FINAL, 
                              replic, which_pops, params, 
                              ld, sitreps, sero, 
                              virus, sgtf, set_id, 
                              opt_mobility, seas_yn, seas_amp, 
                              opt_v2, opt_relu, voct, voci, vacc) {

  postI <- base_run$posteriorsI
  paramsI_base <- base_run$parametersI
  priorsI <- base_run$priorsI
  test_base <- base_run$test
  constants <- base_run$constants
  
  if(vacc_scenario != "trade-off"){ 
  params <- paramsI_base[[1]]
  }
# Counterfactual scenario options
  
  if(vacc_scenario == "no_vaccination"){
    for(i in 1:length(vacc[[1]]$v)){
      vacc[[1]]$v[[i]] <- rep(0,16)
    }
  }
  
  
  if(vacc_scenario == "pfizer-instead"){
    
    # Change vaccine efficacy parameters to Pfizer
    
    params$pop[[1]]$ei_v <- rep(0.85,16)
    ed_oi = 0.9
    params$pop[[1]]$ei2_v <- rep(0.8,16)
    ed_oi2 = 0.81
    params$pop[[1]]$ed_vi <- calc_ve_d(params$pop[[1]]$ei_v,ed_oi)
    params$pop[[1]]$ed_vi2 <- calc_ve_d(params$pop[[1]]$ei2_v,ed_oi2)
  }
  
  if(vacc_scenario == "pfizer-delay"){
    
    # Change vaccine efficacy parameters to Pfizer
    
    params$pop[[1]]$ei_v <- rep(0.85,16)
    ed_oi = 0.9
    params$pop[[1]]$ei2_v <- rep(0.8,16)
    ed_oi2 = 0.81
    params$pop[[1]]$ed_vi <- calc_ve_d(params$pop[[1]]$ei_v,ed_oi)
    params$pop[[1]]$ed_vi2 <- calc_ve_d(params$pop[[1]]$ei2_v,ed_oi2)
    
    # Add delay to vaccination campaign - 2 months later starting mid-April 2021.
    
    for(i in 1:length(vacc[[1]]$v)){
      vacc[[1]]$vt[[i]] <- vacc[[1]]$vt[[i]] + (7*8)
    }
  }
  
  if(vacc_scenario == "az-instead"){
    
    # Change vaccine efficacy parameters to AZ  
    
    params$pop[[1]]$ei_v <- rep(0.75,16)
    params$pop[[1]]$ei2_v <- rep(0.63,16)
    ed_oi <- rep(0.8,16)
    ed_oi2 <- rep(0.65,16)
    params$pop[[1]]$ed_vi <- calc_ve_d(params$pop[[1]]$ei_v,ed_oi)
    params$pop[[1]]$ed_vi2 <- calc_ve_d(params$pop[[1]]$ei2_v,ed_oi2)
    
  }
  
  if(vacc_scenario == "az-delay"){
    
    # Change vaccine efficacy parameters to AZ
    
    params$pop[[1]]$ei_v <- rep(0.75,16)
    params$pop[[1]]$ei2_v <- rep(0.63,16)
    ed_oi <- rep(0.8,16)
    ed_oi2 <- rep(0.65,16)
    params$pop[[1]]$ed_vi <- calc_ve_d(params$pop[[1]]$ei_v,ed_oi)
    params$pop[[1]]$ed_vi2 <- calc_ve_d(params$pop[[1]]$ei2_v,ed_oi2)
    
    
    # Add delay to vaccination campaign - 2 months later starting mid-April 2021
    
    for(i in 1:length(vacc[[1]]$v)){
      vacc[[1]]$vt[[i]] <- vacc[[1]]$vt[[i]] + (7*8)
    }
  }
  

  posteriorsI <- list()
  dynamicsI <- list()
  parametersI <- list()
  
  # Check execution time
   start_time <- Sys.time()
    
    # Loop through regions
    for (pn in which_pops) {
      paramsI = rlang::duplicate(params);
      paramsI$pop = list(rlang::duplicate(params$pop[[pn]]));
      paramsI$travel = matrix(1, nrow = 1, ncol = 1);
      paramsI$schedule = list();
      j = 1;
      
      for (i in seq_along(params$schedule)) {
        if (pn - 1 == params$schedule[[i]]$pops) {
          paramsI$schedule[[j]] = rlang::duplicate(params$schedule[[i]]);
          paramsI$schedule[[j]]$pops = 0;
          j = j + 1;
        }
      }
      
    
      ldI <- rlang::duplicate(ld);
      ldI <- ldI[pid == pn - 1];
      sitrepsI <- rlang::duplicate(sitreps);
      sitrepsI <- sitrepsI[pid == pn - 1];
      seroI <- rlang::duplicate(sero);
      seroI <- seroI[pid == pn - 1];   # sero: all but NHSBT
      virusI <- rlang::duplicate(virus);
      virusI <- virusI[pid == pn - 1 & Data.source %like% "REACT"]; # virus: REACT only
      sgtfI <- copy(sgtf);
      sgtfI <- sgtfI[pid == pn - 1];
      
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
      posteriorsI[[pn]] <- postI[[1]]
      parametersI[[pn]] <- rlang::duplicate(paramsI)
      print(pn)
    }
    
    # Check time
    time2 <- Sys.time()
    print(time2-start_time)

    # Sample dynamics from fit
    
    dynamicsI = list()
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
      paramsI2 = rlang::duplicate(parametersI[[pn]])
      paramsI2$time1 = as.character(ymd(parametersI[[pn]]$time1) + 56);
      test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriorsI[[pn]], 500, seed = 0); # NB 500 is number of samples
      rows = cm_backend_sample_fit_rows(cm_translate_parameters(paramsI2), posteriorsI[[pn]], 500, seed = 0);
      
      test = rbindlist(test)
      test[, population := pn]
      
      # Add dispersion parameters
      disp = posteriorsI[[pn]][rows, .SD, .SDcols = patterns("^disp|v2_conc|v2_disp|v2_sgtf0")]
      disp[, run := .I]
      test = merge(test, disp, by = "run")
      
      dynamicsI[[pn]] = test
      
    }
    
    
    # Concatenate dynamics 
    test <- rbindlist(dynamicsI, fill = TRUE)
    test$population <- as.character(test$population)
    test[, population := "Dominican Republic"]
    
  # Calculate output for baseline and counterfactual (burden averted)
  options(scipen=999) 
  
  counterfactual_fit <- calc_all_averted_burden(test_base, test)
  print(Sys.time() - start_time)
  return(counterfactual_fit)
}

# Function to get tables of burden averted and burden under the counterfactual scenario, baseline scenario (original model fit)
# and from data for comparison

get_burden_tables <- function(d, fit, sim_end){
  end_date <- sim_end
  data <- d
  
  burden_averted_table_6 <- fit %>% 
    filter(ModelType == "Averted" & AgeBand == "All" & (ValueType == "deaths" | ValueType == "hospital_beds" | ValueType == "hosp_ad" | ValueType == "icu_beds" | ValueType == "icu_i")) %>% 
    filter(date  >= as.Date("2021-02-16")) %>% 
    filter(date <= as.Date( "2021-08-16")) %>% 
    group_by(ValueType) %>%
    summarise(total_averted = sum(Value),  total_05 = sum(`Quantile 0.05`), total_25 = sum(`Quantile 0.25`), total_75 = sum(`Quantile 0.75`), total_95 = sum(`Quantile 0.95`)) %>% 
    mutate(run = "Burden Averted")
  
  
  scenario_table_6 <- fit %>% 
    filter(ModelType == "Scenario" & AgeBand == "All" & (ValueType == "deaths" | ValueType == "hospital_beds" | ValueType == "hosp_ad" | ValueType == "icu_beds" | ValueType == "icu_i")) %>% 
    filter(date  >= as.Date("2021-02-16")) %>% 
    filter(date <= as.Date( "2021-08-16")) %>% 
    group_by(ValueType) %>%
    summarise(total_averted = sum(Value),  total_05 = sum(`Quantile 0.05`), total_25 = sum(`Quantile 0.25`), total_75 = sum(`Quantile 0.75`), total_95 = sum(`Quantile 0.95`)) %>% 
    mutate(run = "Scenario")
  
  
  base_table_6 <- fit %>% 
    filter(ModelType == "Baseline" & AgeBand == "All" & (ValueType == "deaths" | ValueType == "hospital_beds" | ValueType == "hosp_ad" | ValueType == "icu_beds" | ValueType == "icu_i")) %>% 
    filter(date  >= as.Date("2021-02-16")) %>% 
    filter(date <= as.Date( "2021-08-16")) %>% 
    group_by(ValueType) %>%
    summarise(total_averted = sum(Value),  total_05 = sum(`Quantile 0.05`), total_25 = sum(`Quantile 0.25`), total_75 = sum(`Quantile 0.75`), total_95 = sum(`Quantile 0.95`)) %>% 
    mutate(run = "Baseline")
  
  
  burden_averted_table_total <- fit %>% 
    filter(ModelType == "Averted" & AgeBand == "All" & (ValueType == "deaths" | ValueType == "hospital_beds" | ValueType == "hosp_ad" | ValueType == "icu_beds" | ValueType == "icu_i")) %>% 
    filter(date  >= as.Date("2021-02-16")) %>% 
    filter(date <= as.Date(end_date)) %>% 
    group_by(ValueType) %>%
    summarise(total_averted = sum(Value),  total_05 = sum(`Quantile 0.05`), total_25 = sum(`Quantile 0.25`), total_75 = sum(`Quantile 0.75`), total_95 = sum(`Quantile 0.95`)) %>% 
    mutate(run = "Burden Averted")
  
  
  scenario_table_total <- fit %>% 
    filter(ModelType == "Scenario" & AgeBand == "All" & (ValueType == "deaths" | ValueType == "hospital_beds" | ValueType == "hosp_ad" | ValueType == "icu_beds" | ValueType == "icu_i")) %>% 
    filter(date  >= as.Date("2021-02-16")) %>% 
    filter(date <= as.Date(end_date)) %>% 
    group_by(ValueType) %>%
    summarise(total_averted = sum(Value),  total_05 = sum(`Quantile 0.05`), total_25 = sum(`Quantile 0.25`), total_75 = sum(`Quantile 0.75`), total_95 = sum(`Quantile 0.95`)) %>% 
    mutate(run = "Scenario")
  
  
  base_table_total <- fit %>% 
    filter(ModelType == "Baseline" & AgeBand == "All" & (ValueType == "deaths" | ValueType == "hospital_beds" | ValueType == "hosp_ad" | ValueType == "icu_beds" | ValueType == "icu_i")) %>% 
    filter(date  >= as.Date("2021-02-16")) %>% 
    filter(date <= as.Date(end_date)) %>% 
    group_by(ValueType) %>%
    summarise(total_averted = sum(Value),  total_05 = sum(`Quantile 0.05`), total_25 = sum(`Quantile 0.25`), total_75 = sum(`Quantile 0.75`), total_95 = sum(`Quantile 0.95`)) %>% 
    mutate(run = "Baseline")
  
  
  data_table <- data %>% 
    filter(ValueType == "Deaths" | ValueType == "Hospital beds\noccupied" | ValueType == "ICU beds\noccupied") %>% 
    filter(d >= as.Date("2021-02-16")) %>% 
    filter(d <= as.Date(end_date)) %>% 
    mutate(Year = year(d)) %>% 
    group_by(ValueType) %>% 
    summarise(total = sum(y)) %>% 
    mutate(Run = "Data") %>% 
    mutate(Year = "All")
  
  return(list(burden_averted_table_6 = burden_averted_table_6, scenario_table_6 = scenario_table_6, base_table_6 = base_table_6, burden_averted_table_total = burden_averted_table_total, scenario_table_total = scenario_table_total, base_table_total = base_table_total))

  }



