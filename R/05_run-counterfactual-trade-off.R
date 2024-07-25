# -------------------------------------------------------------------------------------------------------
# Run counterfactual analysis scenarios looking at the trade-off between vaccination coverage
# and population mobility
# Author: Emilie Finch
# -------------------------------------------------------------------------------------------------------

# Set end date for counterfactual simulations 
## Here comparing burden in 6 months following start of vaccination campaign

date_end = "2021-08-16" 
today <- Sys.Date()

# Get list of mobility changes  -------------------------------------------------------------------------
## Work and other contacts reduced from -40%:40% in 10% intervals from 16th Feb 2021 until 16th August 2021

list_schedules <- list()
intervals <- seq(-0.4,0.4, by = 0.1)

for(k in intervals){
  new_schedule <- schedule[[1]]
  # For each day after vacc campaign start and before simulation end
  for (j in 412:(as.Date(date_end) - as.Date("2020-01-01"))) { 
    # Only change work and other contacts (not home or school)
    for(q in 2:5){ 
      new_schedule$values[[j]][q] <- max(schedule[[1]]$values[[j]][q]-k,0)
    }
  }
  list_schedules[[which(k == intervals)]] <- new_schedule
}

# Get list of new vaccination schedules -----------------------------------------------------------------
## Vaccination coverage from 0 - 100% by date_end

original_vacc <- ox_vacc %>% 
  filter(date <= date_end) %>% 
  select(people_new_fully_vaccinated_approx) %>% 
  pull()

list_vacc_schedules <- list()
proportions <- seq(0, 1, by = 0.1)

plots <- list()
for(m in 1:length(proportions)){
  
  # Allocate extra vaccines by multiplying daily Oxford vaccine data time series by the proportion change
  # between the final reported proportion vaccinated and the simulation final proportion
  # and then allocating doses between eligible age groups using the same process as when building the vaccination data for fitting
  
  new_vacc <- allocate_extra_vaccines(original_vacc, eligible_ages, age_pops, date_end, proportions[m])
  allocation_plots <- plot_vacc_allocation(new_vacc, age_pops, proportions[m])
  list_vacc_schedules[[m]] <- list(vt = ox_vacc$date[1:(as.Date(date_end) - min(ox_vacc$date))], v = new_vacc$vaccs_age)
  plots[[m]] <- allocation_plots$perc_plot
}

# Run counterfactual for each combination of mobility change and vaccination coverage ------------------

burden_out <- NULL

list_counterfactual_out <- list()

for(o in 1:length(list_schedules)){ 
  params$schedule <- list()
  params$schedule[[1]] <- list_schedules[[o]]
  
  for(p in 1:length(list_vacc_schedules)){ 
    
    new_vacc <- vacc
    new_vacc[[1]] <- list_vacc_schedules[[p]]
    
    cat(paste0("Running scenarios: contact option ", o, "/", length(list_schedules), " ,  coverage option ", p, "/", length(list_vacc_schedules)))  
    scenario_name = paste0("scenario_contacts-", intervals[o], "_vacc-coverage-", proportions[p])
    
    
  # Run counterfactual scenario with new vaccine coverage and contact schedule 
    
   counterfactual_fit <- run_counterfactual(base_run = mod_run, sim_end = "2021-12-15", vacc_scenario = "trade-off",
                                                 REP_START, REP_END, BURN_IN, 
                                                 BURN_IN_FINAL, ITER, ITER_FINAL, 
                                                 replic, which_pops,
                                                 params, ld, sitreps, sero, 
                                                 virus, sgtf, set_id, opt_mobility, 
                                                 seas_yn, seas_amp, opt_v2, 
                                                 opt_relu, voct, voci, vacc =new_vacc)
    
    options(scipen=999) 
    comp_outs <- get_burden_tables(counterfactual_fit, d = mod_fit[[1]], sim_end = date_end) 
    comp_outs <- c(comp_outs, label = scenario_name)
    
   # Total
    
    add_total <- data.frame(contact = intervals[o], 
                      vaccination = proportions[p],
                      ValueType = comp_outs$burden_averted_table$ValueType[comp_outs$burden_averted_table$horizon == "Total"],
                      value = comp_outs$burden_averted_table$total_averted[comp_outs$burden_averted_table$horizon == "Total"],
                      total_05 = comp_outs$burden_averted_table$total_05[comp_outs$burden_averted_table$horizon == "Total"],
                      total_25 = comp_outs$burden_averted_table$total_25[comp_outs$burden_averted_table$horizon == "Total"],
                      total_75 = comp_outs$burden_averted_table$total_75[comp_outs$burden_averted_table$horizon == "Total"],
                      total_95 = comp_outs$burden_averted_table$total_95[comp_outs$burden_averted_table$horizon == "Total"])
    
    burden_out <- bind_rows(burden_out, add_total)
    
    qsave(list(burden_table = burden_out), file = here("output", paste0("counterfactual-trade-off_", today, ".qs")))
    
  }
}


qsave(list(burden_table = burden_out), file = here("output", paste0("counterfactual-trade-off_", today, ".qs")))
