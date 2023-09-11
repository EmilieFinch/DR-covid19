
# Function to generate age-stratified vaccination data from daily total vaccinations, eligible age groups and population sizes

# vaccs (numeric): number of people fully vaccinated every day
# eligible ages (list): which age groups are eligible for vaccination at each time point (1 = eligible, 0 = not eligible)
# age_pops (numeric): population size of each age group

allocate_vaccines <- function(vaccs, eligible_ages, age_pops) { 

vaccs_age <- NULL
groups_vaccinated <- NULL
tracking <- rep(0,16)

for (i in 1:length(vaccs)) {
groups_vaccinated[[i]] <- eligible_ages[[i]]

for (age_group in 1:length(age_pops)) {
  if (tracking[age_group] >= age_pops[age_group]) {
    groups_vaccinated[[i]][age_group] <- 0 
    # If eligible age group has already been totally vaccinated then they are not included in groups_vaccinated
  }
}

num_groups <- sum(groups_vaccinated[[i]]) # Number of age groups vaccinated at this time point

if (num_groups > 0) { 
vaccs_age[[i]] <- groups_vaccinated[[i]]*(vaccs[i]/num_groups)
} else {
  print(paste0("All eligible groups fully vaccinated at time step ", i ))
# If age groups eligible for vaccination are already totally vaccinated, spread available vaccines at that time point 
# between other age groups > 20. 
# This simulates some individuals in other age groups being vaccinated due to HCW status or other conditions.
  new_num <- 12 - sum(eligible_ages[[i]])
  vaccs_age[[i]] <- vaccs[i]/new_num * (1-eligible_ages[[i]])*rep(c(0,1), times = c(4,12))
}

int_tracking <- tracking + vaccs_age[[i]]

# If an age group will be more than fully vaccinated in this time step, reallocate the vaccines again between all adult age groups
surplus_doses <- 0
ages_for_surplus <- groups_vaccinated[[i]]

for (age_group in 1:length(age_pops)) {
  
  # If total vaccinated at the end of this iteration in an age group is greater than the population size, only vaccinate with enough
  # doses to match population size and reallocate other doses among other adult age groups not yet totally vaccinated
  
  if (int_tracking[age_group] > age_pops[age_group]) { 
    print(paste0("Surplus doses in age group ", age_group, " at time step ", i , " reallocating to other groups"))
    ages_for_surplus[age_group] <- 0
    surplus_doses <- surplus_doses + (int_tracking[age_group] - age_pops[age_group])
    vaccs_age[[i]][age_group] <- vaccs_age[[i]][age_group] - surplus_doses
  }
}

# Allocate to surplus doses to rest of adult population not already targeted

if (surplus_doses > 0 & sum(eligible_ages[[i]]) < 12) {
  new_num <- 12 - sum(eligible_ages[[i]])
  vaccs_age[[i]] <- vaccs_age[[i]] + surplus_doses/new_num * (1-eligible_ages[[i]])*rep(c(0,1), times = c(4,12))
  
  # If all adult population already targeted allocate surplus doses between these again
  
} else if(surplus_doses > 0 & sum(eligible_ages[[i]]) == 12){
  new_num <- sum(ages_for_surplus)
  vaccs_age[[i]] <- vaccs_age[[i]] + surplus_doses/new_num * ages_for_surplus
  
}

tracking <- tracking + vaccs_age[[i]]
}
return(list(vaccs_age = vaccs_age, groups_vaccinated = groups_vaccinated))
}

# Function to generate age-stratified vaccination data for counterfactual scenarios, to achieve a specified vaccination coverage 

# vaccs (numeric): number of people fully vaccinated every day
# eligible ages (list): which age groups are eligible for vaccination at each time point (1 = eligible, 0 = not eligible)
# age_pops (numeric): population size of each age group
# date_end

allocate_extra_vaccines <- function(ox_vaccs, eligible_ages, age_pops, date_end, sim_coverage) {
 
   ox_vaccs <- ox_vaccs[1:(as.Date(date_end) - as.Date("2021-02-15"))]
  rep_coverage <- sum(ox_vaccs) / sum(age_pops) # Reported final coverage by end of simulation period
  percentage_change <- sim_coverage / rep_coverage
  vaccs_age <- NULL
  groups_vaccinated <- NULL
  young_vaccinated <- NULL
  tracking <- rep(0,16)
  
  for(i in 1:43){
    vaccs_age[[i]] <- rep(0,16)
  }
  
  # From first day of data for individuals fully vaccinated
  for(i in 44:length(ox_vaccs)) { 
    groups_vaccinated[[i]] <- eligible_ages[[i]]
    young_vaccinated[[i]] <- rep(c(1,0), times = c(4,12))
    
    for(age_group in 1:length(age_pops)){
      if(tracking[age_group] >= age_pops[age_group]){
        # If eligible age group has already been totally vaccinated then they are not included in groups_vaccinated
        groups_vaccinated[[i]][age_group] <- 0 
        young_vaccinated[[i]][age_group] <- 0
      }
    }
    
    # Number of eligible age groups vaccinated at this time point
    num_groups <- sum(groups_vaccinated[[i]]) 

      if(num_groups > 0) { 
        
      vaccs_age[[i]] <-  pmax(groups_vaccinated[[i]]*(ox_vaccs[i]*percentage_change/num_groups),0)
      
    } else if (num_groups == 0 & i < 83) {
      
      print(paste0("All eligible groups fully vaccinated at time step ", i, ", distributing vaccines amongst other ages > 20"))
      # If age groups eligible for vaccination are already totally vaccinated, spread available vaccines at that time point between other age groups > 20. 
      # This simulates some individuals in other age groups being vaccinated due to HCW status or other conditions.
      new_num <- 12 - sum(eligible_ages[[i]])
      vaccs_age[[i]] <- pmax((ox_vaccs[i]*percentage_change/new_num) * (1-eligible_ages[[i]])*rep(c(0,1), times = c(4,12)),0)
      
    } else if(num_groups == 0 & i > 83 & sum(tracking >= age_pops) < 16) {
      
      print(paste0("All adults fully vaccinated at time step ", i, ", distributing vaccines amongst other ages < 20"))
      num_young <- sum(young_vaccinated[[i]])
      vaccs_age[[i]] <- pmax((ox_vaccs[i]*percentage_change/num_young) * young_vaccinated[[i]],0)
      groups_vaccinated[[i]] <- young_vaccinated[[i]]
      
    } else if(sum(tracking >= age_pops) == 16){
      
      print(paste0("Entire population vaccinated at time step ", i, ", no further vaccination."))
      vaccs_age[[i]] <- rep(0,16)
    }
    
    int_tracking <- tracking + vaccs_age[[i]]
    
    # If an age group will be more than fully vaccinated in this time step, reallocate the vaccines again
    surplus_doses <- 0
    ages_for_surplus <- groups_vaccinated[[i]]
    young_for_surplus <- rep(c(1,0), times = c(4,12))
    
    for(age_group in 1:length(age_pops)) {
      # If total vaccinated at the end of this iteration in an age group is greater than the population size, only vaccinate with enough doses to match population size
      # And reallocate other doses among other adult age groups not yet totally vaccinated
      
      if(int_tracking[age_group] > age_pops[age_group]) { 
        print(paste0("Surplus doses in age group ", age_group, " at time step ", i , " reallocating to other groups"))
        ages_for_surplus[age_group] <- 0
        young_for_surplus[age_group] <- 0
        
        surplus_doses <- surplus_doses + (int_tracking[age_group] - age_pops[age_group])
        vaccs_age[[i]][age_group] <- vaccs_age[[i]][age_group] - (int_tracking[age_group] - age_pops[age_group])
      }
    }
    
    # Allocate to surplus doses to rest of adult population not already targeted
    if(surplus_doses > 0 & sum(eligible_ages[[i]]) < 12) {
      new_num <- 12 - sum(eligible_ages[[i]])
      vaccs_age[[i]] <- vaccs_age[[i]] + surplus_doses/new_num * (1-eligible_ages[[i]])*rep(c(0,1), times = c(4,12))
      print(paste0("Surplus doses reallocated to other ages > 20"))
      
      # If all adult population already targeted allocate surplus doses between these again
    } else if(surplus_doses > 0 & sum(eligible_ages[[i]]) == 12 & sum(ages_for_surplus > 1)){
      new_num <- sum(ages_for_surplus)
      vaccs_age[[i]] <- vaccs_age[[i]] + surplus_doses/new_num * ages_for_surplus
      print(paste0("Surplus doses reallocated to eligible age groups"))
      
    } else if(surplus_doses > 0 & sum(eligible_ages[[i]]) == 12 & sum(ages_for_surplus <=1) & sum(young_for_surplus > 0)){
      new_num <- sum(young_for_surplus)
      vaccs_age[[i]] <- vaccs_age[[i]] + surplus_doses/new_num * young_for_surplus
      print(paste0("Surplus doses reallocated to < 20 year olds"))
    } else if(surplus_doses >0 & sum(eligible_ages[[i]] == 12 & sum(ages_for_surplus == 0))){
      print(paste0("Entire population vaccinated in this time step, surplus doses not reallocated"))
      
    } else {}
    
        tracking <- tracking + vaccs_age[[i]]
  }
  
  return(list(vaccs_age = vaccs_age, groups_vaccinated = groups_vaccinated))
  
}

plot_vacc_allocation <- function(vaccs_out, age_pops, final_coverage) {
  
  # Check output
  total_vaccs <- rep(0,16)
  prop_vaccinated <- NULL
  df_vaccinated <- data.frame(time_step = numeric(),
                              age_group = numeric(),
                              percentage = numeric(),
                              new_vaccs = numeric(),
                              cumulative_vaccs = numeric())
  
  df <- c()
  for(v in 1:length(vaccs_out$vaccs_age)) {
    total_vaccs <- total_vaccs + vaccs_out$vaccs_age[[v]]
    prop_vaccinated[[v]] <- total_vaccs/age_pops*100
    for(age in 1:length(age_pops)){ 
      df$time_step <- as.numeric(v)
      df$age_group <- as.numeric(age)
      df$percentage <- prop_vaccinated[[v]][age]
      df$new_vaccs <- vaccs_out$vaccs_age[[v]][age]
      df$cumulative_vaccs <- total_vaccs[age]
      df_vaccinated <- rbind(df_vaccinated, df)
    }
  }
  
  df_vaccinated <- as.data.frame(df_vaccinated)
  
  # Plot percentage vaccinated by age group to check
  age_labs <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+")
  df_vaccinated$age_group <- factor(df_vaccinated$age_group, labels = age_labs)
  
  perc_plot <- ggplot(df_vaccinated) +
    geom_line(aes(x = as.Date("2021-02-15") + time_step, y = percentage, group = factor(age_group), col = factor(age_group))) +
    geom_hline(aes(yintercept = final_coverage*100), lty = "dashed") + 
    scale_y_continuous(breaks = seq(0,110, by = 10)) +
    labs(x = "Date", y = "Percentage of age group fully vaccinated", col = "Age Group")
  
  # Plot new vaccinated by age group to check 
  
  new_vaccinated_plot <- ggplot(df_vaccinated) +
    geom_line(aes(x = as.Date("2021-02-15") + time_step, y = new_vaccs, group = factor(age_group), col = factor(age_group))) +
    labs(x = "Date", y = "New people fully vaccinated", col = "Age Group")
  
  
  # Plot cumulative vaccinated by age group to check
  
  cumulative_vaccinated_plot <- ggplot(df_vaccinated) +
    geom_line(aes(x = as.Date("2021-02-15") + time_step, y = cumulative_vaccs, group = factor(age_group), col = factor(age_group))) +
    labs(x = "Date", y = "Cumulative people fully vaccinated", col = "Age Group")
  
  
  # Calculate total vaccines to check
  
  total_allocated <- sum(df_vaccinated %>% 
        group_by(age_group) %>% 
        summarise(total = max(cumulative_vaccs)) %>% 
        pull(total))
  
  return(list(perc_plot = perc_plot, new_vaccinated_plot = new_vaccinated_plot, cumulative_vaccinated_plot = cumulative_vaccinated_plot, total_allocated = total_allocated))
}
