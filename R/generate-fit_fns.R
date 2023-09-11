# Functions to generate model fits from posteriors

# Function to wrangle dynamics and add dispersion from observation process

wrangle_dynamics <- function(summ, varname, outname, cyear, cmonth, cday, ymd_from, ymd_to, size)
{
  quants <- seq(0.05, 0.95, by = 0.05)
  
  rq <- data.table(run = 1:summ[, max(run)], q = runif(summ[, max(run)]))
  summ <- merge(summ, rq, by = "run")
  if (is.character(size)) {
    summ[, rv := qnbinom(q, size = 1.0 / get(size)^2, mu = get(varname))]
  } else if (size == 0) { # poisson
    summ[, rv := qpois(q, lambda = get(varname))]
  } else if (size == -1) { # no variation
    summ[, rv := get(varname)]
  } else if (size > 0) { # neg binom
    summ[, rv := qnbinom(q, size = size, mu = get(varname))]
  } else {
    stop("Size must be -1, 0, or a positive number")
  }
  qsumm <- summ[, as.list(quantile(rv, quants, na.rm = TRUE)), by = .(t, population, age_group)]
    qsumm[, date := t + ymd("2020-01-01")]
  qsumm <- rbind(qsumm,
                qsumm[, lapply(.SD, sum), .SDcols = `5%`:`95%`, by = .(t, population, date, age_group = rep("All", nrow(qsumm)))],
                use.names = TRUE, fill = TRUE)
  
  data_out <- qsumm[date %between% c(ymd_from, ymd_to),
                   .(Group = "LSHTM", Model = "Transmission", Scenario = "MTP", ModelType = "Multiple", Version = 2,
                     `Creation Day` = cday, `Creation Month` = cmonth, `Creation Year` = cyear,
                     `Day of Value` = day(date), `Month of Value` = month(date), `Year of Value` = year(date),
                     AgeBand = age_group, Geography = population, ValueType = outname, Value = `50%`,
                     `5%`, `10%`, `15%`, `20%`, `25%`, `30%`, `35%`, `40%`, `45%`, `50%`, `55%`, `60%`, `65%`, `70%`, `75%`, `80%`, `85%`, `90%`, `95%`)]
  
  names(data_out)[names(data_out) %like% "\\%$"] = paste("Quantile", quants);
  data_out
}

# Function to combine wrangled dynamics for all variables

combine_dynamics <- function(test0, cyear, cmonth, cday, ymd_from, ymd_to)
{
  cat("Restricting time span...\n")
  t0 <- as.numeric(ymd(ymd_from) - ymd("2020-01-01"))
  t1 <- as.numeric(ymd(ymd_to) - ymd("2020-01-01"))
  test <- test0[t %between% c(t0, t1)]
  test <- test0
  
  cat("Adding age groups...\n");
  test[group %between% c(1, 1), age_group := "0-4"]
  test[group %between% c(2, 3), age_group := "5-14"]
  test[group %between% c(4, 5), age_group := "15-24"]
  test[group %between% c(6, 9), age_group := "25-44"]
  test[group %between% c(10, 13), age_group := "45-64"]
  test[group %between% c(14, 15), age_group := "65-74"]
  test[group %between% c(16, 16), age_group := "75+"]
  test[, age_group := factor(age_group, levels = unique(age_group))]
  
  cat("Summarizing variables...\n");
  summ <- test[, .(death_o = sum(death_o + death2_o), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
  summ3 <- test[, .(icu_p = sum(icu_p + icu2_p), disp = mean(disp_icu_prev)), by = .(run, t, population, age_group)]
  summ4 <- test[, .(bed_p = sum(pmax(0, hosp_p + hosp2_p - hosp_undetected_p - hosp_undetected2_p)), disp = mean(disp_hosp_prev)), by = .(run, t, population, age_group)]
  summ34 <- test[, .(admissions = sum(hosp_undetected_o + hosp_undetected2_o), disp = disp_hosp_inc), by = .(run, t, population, age_group)]
  summ_inf_i <- test[, .(infections_i = sum(pcr_positive_i)), by = .(run, t, population, age_group)]
  summ_inf_p <- test[, .(infections_p = sum(pcr_positive_p)), by = .(run, t, population, age_group)];
  summ_sero_p <- test[, .(sero_p = sum(lfia_positive_p)), by = .(run, t, population, age_group)];
  summ_rt <- test[group == 1, .(age_group = "All", Rt = obs0), by = .(run, t, population)]
  summ_r0 <- test[group == 4, .(age_group = "All", R0 = obs0), by = .(run, t, population)]
  summ_rt1 <- test[group == 11, .(age_group = "All", Rt_1 = obs0), by = .(run, t, population)]
  summ_rt2 <- test[group == 12, .(age_group = "All", Rt_2 = obs0), by = .(run, t, population)]
  
  
  cat("Running quantiles...\n");
  w <- rbind(
    wrangle_dynamics(summ, "death_o", "type28_death_inc_line", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
    wrangle_dynamics(summ3, "icu_p", "icu_prev", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
    wrangle_dynamics(summ4, "bed_p", "hospital_prev", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
    wrangle_dynamics(summ34, "admissions", "hospital_inc", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
    wrangle_dynamics(summ_inf_p, "infections_p", "prevalence_mtp", cyear, cmonth, cday, ymd_from, ymd_to, -1),
    wrangle_dynamics(summ_inf_i, "infections_i", "infections_inc", cyear, cmonth, cday, ymd_from, ymd_to, -1),
    wrangle_dynamics(summ_sero_p, "sero_p", "sero_prev", cyear, cmonth, cday, ymd_from, ymd_to, -1),
    wrangle_dynamics(summ_rt, "Rt", "Rt", cyear, cmonth, cday, ymd_from, ymd_to, -1),
    wrangle_dynamics(summ_r0, "R0", "R0", cyear, cmonth, cday, ymd_from, ymd_to, -1),
    wrangle_dynamics(summ_rt1, "Rt_1", "Rt_1", cyear, cmonth, cday, ymd_from, ymd_to, -1),
    wrangle_dynamics(summ_rt2, "Rt_2", "Rt_2", cyear, cmonth, cday, ymd_from, ymd_to, -1)
    
  )
  return (w)
}

make_data <- function(ld, sitreps, virus, sero)
{
  rbind(
    ld[, .(ValueType = "type28_death_inc_line", Geography = name,
           dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = N, ymax = NA)],
    sitreps[, .(ValueType = "icu_prev", Geography = name,
                dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_in_itu, ymax = NA)],
    sitreps[, .(ValueType = "hospital_prev", Geography = name,
                dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_in_all_beds, ymax = NA)],
    sitreps[, .(ValueType = "hospital_inc", Geography = name,
                dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_admitted_diagnosed, ymax = NA)],
    virus[Data.source %like% "REACT", .(ValueType = "prevalence_mtp", Geography = NHS.region,
                                        dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
                                        ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))],
    sero[, .(ValueType = "sero_prev", Geography = NHS.region,
             dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
             ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))]
  )
}

# Wrapper function to generate fit from MCMC model trajectory output (wrangle and combine dynamics)

gen_fit <- function(test, parametersI, ld, sitreps, virus, sero, populations)
{
    test <- copy(test)
    sero <- copy(sero)
    sero <- sero[NHS.region %in% populations]
    virus <- copy(virus)
    virus <- virus[NHS.region %in% populations]
    ld <- copy(ld)
    ld <- ld[name %in% populations]
    sitreps <- copy(sitreps)
    sitreps <- sitreps[name %in% populations]
    
    # Calculate total population
    popsize <- NULL
    for (i in seq_along(parametersI)) {
        if (!is.null(parametersI[[i]])) {
            popsize = rbind(popsize,
                data.table(Geography = parametersI[[i]]$pop[[1]]$name, population_size = sum(parametersI[[i]]$pop[[1]]$size))
            )
        }
    }
    
    # Create formatted output
    
    today <- Sys.Date()
    year <- year(today)
    month <- month(today)
    day <- day(today)
    
    output <- combine_dynamics(test, year, month, day, "2020-01-01", as.character(ymd("2020-01-01") + max(test$t)))
    output[, d := make_date(`Year of Value`, `Month of Value`, `Day of Value`)]
    output <- merge(output, popsize, by = "Geography")
    
    adj_output <- function(output, val_type, div, pop = 0) {
        output[ValueType == val_type, Value := Value / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.05` := `Quantile 0.05` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.1`  := `Quantile 0.1`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.15` := `Quantile 0.15` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.2`  := `Quantile 0.2`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.25` := `Quantile 0.25` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.3`  := `Quantile 0.3`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.35` := `Quantile 0.35` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.4`  := `Quantile 0.4`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.45` := `Quantile 0.45` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.5`  := `Quantile 0.5`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.55` := `Quantile 0.55` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.6`  := `Quantile 0.6`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.65` := `Quantile 0.65` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.7`  := `Quantile 0.7`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.75` := `Quantile 0.75` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.8`  := `Quantile 0.8`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.85` := `Quantile 0.85` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.9`  := `Quantile 0.9`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.95` := `Quantile 0.95` / (div + population_size * pop)]
    }
    
    adj_output(output, "hospital_inc", 1)
    adj_output(output, "hospital_prev", 1)
    adj_output(output, "icu_prev", 1)
    adj_output(output, "prevalence_mtp", 0, 0.01)
    adj_output(output, "sero_prev", 0, 0.01)
    adj_output(output, "type28_death_inc_line", 1)
    
    # Make data to output
    
    data <- make_data(ld, sitreps, virus, sero)
    data <- merge(data, popsize, by = "Geography")
    
    adj_data <- function(data, val_type, div, pop = 0) {
        data[ValueType == val_type, ymin := ymin / (div + population_size * pop)]
        data[ValueType == val_type, y    := y    / (div + population_size * pop)]
        data[ValueType == val_type, ymax := ymax / (div + population_size * pop)]
    }
    
    adj_data(data, "hospital_inc", 1)
    adj_data(data, "hospital_prev", 1)
    adj_data(data, "icu_prev", 1)
    adj_data(data, "prevalence_mtp", 0.01)
    adj_data(data, "sero_prev", 0.01)
    adj_data(data, "type28_death_inc_line", 1)
    
    output[ValueType == "hospital_inc", ValueType := "Hospital\nadmissions"]
    output[ValueType == "hospital_prev", ValueType := "Hospital beds\noccupied"]
    output[ValueType == "icu_prev", ValueType := "ICU beds\noccupied"]
    output[ValueType == "infections_inc", ValueType := "Infection\nincidence"]
    output[ValueType == "prevalence_mtp", ValueType := "PCR\nprevalence (%)"]
    output[ValueType == "sero_prev", ValueType := "Seroprevalence\n(%)"]
    output[ValueType == "type28_death_inc_line", ValueType := "Deaths"]

    data[ValueType == "hospital_inc", ValueType := "Hospital\nadmissions"]
    data[ValueType == "hospital_prev", ValueType := "Hospital beds\noccupied"]
    data[ValueType == "icu_prev", ValueType := "ICU beds\noccupied"]
    data[ValueType == "infections_inc", ValueType := "Infection\nincidence"]
    data[ValueType == "prevalence_mtp", ValueType := "PCR\nprevalence (%)"]
    data[ValueType == "sero_prev", ValueType := "Seroprevalence\n(%)"]
    data[ValueType == "type28_death_inc_line", ValueType := "Deaths"]
    
    return (list(data, output))
}

