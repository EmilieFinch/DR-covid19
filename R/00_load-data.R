# -------------------------------------------------------------------------------
# Load COVID-19 data for the Dominican Republic
# Author: Emilie Finch
# -------------------------------------------------------------------------------

# This script loads the data and the packages needed for subsequent analyses
# This includes:
# - Surveillance data (number of COVID-19 deaths until 15th December 2021)
# - Hospital and ICU bed occupancy data (number of COVID-19 deaths until 15th December 2021)
# - Serological data (found in Nilles et al, 2022 - https://doi.org/10.1016/j.lana.2022.100390)
# - Vaccination data (downloaded from Oxford's Our World in Data: https://github.com/owid/covid-19-data)
# - Mobility data (downloaded from: https://www.google.com/covid19/mobility/)

library(here)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(lubridate)
library(zoo)
library(qs)
library(ggplot2)
library(cowplot)
library(readxl)
library(sn)
library(mgcv)
library(binom)
library(scales)

source(here("R", "allocate-vaccines_fns.R"))
source(here("R", "utils.R"))

ld <- fread(here("data", "surveillance-deaths_dom-rep.csv")) %>% 
  mutate(date = as.Date(date, format = "%d/%m/%Y"))
sitreps <- fread(here("data", "hospital-icu-beds_dom-rep.csv")) %>% 
  mutate(date = as.Date(date, format = "%d/%m/%Y"))
sero <- fread(here("data", "serology_dom-rep.csv")) %>% 
  mutate(Start.date = as.Date(Start.date, format = "%d/%m/%Y"),
         End.date = as.Date(End.date, format = "%d/%m/%Y"))

ld$date <- as.Date(ld$date, format = "%d/%m/%Y")
sitreps$date <- as.Date(sitreps$date, format = "%d/%m/%Y")

# Add empty virus and STGF data frames (as not available for the DR)

virus <- data.table(NHS.region = NA, Start.date = NA, End.date = NA, Central.estimate = NA, 
                    Lower.bound = NA, Upper.bound = NA, Test = NA, Median.mean = NA, Min.age = NA, 
                    Max.age = NA, Data.source = NA, N.tests = NA, xi = NA, omega = NA, alpha = NA, pid = NA)
sgtf <- data.table(date = NA, nhs_name = NA, sgtf = NA, other = NA, pid = NA)

# Vaccination data

ox_vacc <- read.csv(here("data", "vaccinations_dom-rep.csv"))

# Replace NAs in cumulative totals of total vaccinations, people vaccinated and people fully vaccinated with most recent figure above

ox_vacc <- ox_vacc %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>% 
  select(location, iso_code, date, people_fully_vaccinated) %>% 
  mutate(people_fully_vaccinated = na_if(people_fully_vaccinated, 0)) %>%
  mutate(people_fully_vaccinated_approx = round(approx(people_fully_vaccinated, xout = row_number())$y)) %>% 
  fill(people_fully_vaccinated_approx, .direction = "down") %>% 
  mutate(people_new_fully_vaccinated_approx = people_fully_vaccinated_approx - lag(people_fully_vaccinated_approx)) 
  
ox_vacc$people_new_fully_vaccinated_approx[is.na(ox_vacc$people_new_fully_vaccinated_approx)] <- 0

vaccs_daily <- ox_vacc$people_new_fully_vaccinated_approx

# List of vectors of eligible ages for vaccination at each time step (from VacunateRD: https://vacunate.gob.do/)

eligible_ages <- NULL
for(i in 1:length(vaccs_daily)){
  if(i <= 67){
    eligible_ages[[i]] <- rep(c(0,1), times = c(14,2))
  } else if(i <= 76){
    eligible_ages[[i]] <- rep(c(0,1), times = c(12,2))
  } else if(i <= 83){
    eligible_ages[[i]] <- rep(c(0,1), times = c(10,6))
  }  else {
    eligible_ages[[i]] <- rep(c(0,1), times = c(4,12))
  }
}   

# Population sizes in 5-year age bands
age_pops <- c(1002829, 996519, 977245, 958630, 939563, 914754, 825573, 738809, 668263, 606226, 541435, 472900, 388564, 293575, 204150, 318869)

# Run vaccine allocation
vaccs_out <- allocate_vaccines(vaccs_daily, eligible_ages, age_pops)

# Save final vaccine df
vacc <- list(list(vt = unique(ox_vacc$date), v = vaccs_out$vaccs_age))

rm(vaccs_out)

# Mobility data

googmo <- read.csv(here("data", "google-mobility_dom-rep.csv"))

googmo <- googmo %>% 
  pivot_longer(cols = -c(1:9), names_to = "GPL_TYPE", values_to = "value") %>% 
  mutate(GPL_TYPE = str_replace(GPL_TYPE, "_percent_change_from_baseline", "")) %>% 
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>% 
  mutate(value = 1 + value/100) %>% 
  mutate(t = as.numeric(date - ymd("2020-01-01"))) %>% 
  group_by(GPL_TYPE) %>% 
  # If NAs then fill these here within time series using na.approx by GPL_TYPE
  mutate(value = na.approx(value, na.rm = FALSE)) %>% 
  mutate(smooth_value = rollmean(value, 7, fill = NA)) %>% 
  ungroup() %>% 
  filter(!is.na(smooth_value))

# Reformat google mobility data

googmo <- googmo %>% 
  select(t, GPL_TYPE, smooth_value) %>% 
  pivot_wider(id_cols = t, names_from = GPL_TYPE, values_from = smooth_value) %>% 
  rename(grocpharm = grocery_and_pharmacy, retrec = retail_and_recreation, transit = transit_stations, workplace = workplaces)

# Fill in from Jan 1 2020

googmo <- rbind(googmo, 
           data.frame(t = 0:47, grocpharm = 1, parks = 1, residential = 1, retrec = 1, transit = 1, workplace = 1)) %>% 
  arrange(t)
setDT(googmo)

# Add schools opening and closing

school_close <-  c("2020-03-19", "2021-12-20"); # Can make these lists
school_reopen <- c("2021-09-20", "2022-01-10"); # Schools reopen with start of academic year in 2021.

school_c <- as.numeric(ymd(school_close) - ymd("2020-01-01"))
school_r <- as.numeric(ymd(school_reopen) - ymd("2020-01-01"))

days <- googmo$t

n_closures <- rowSums(matrix(days, ncol = length(school_c), nrow = length(days), byrow = FALSE) >= 
                       matrix(school_c, ncol = length(school_c), nrow = length(days), byrow = TRUE));
n_reopenings <- rowSums(matrix(days, ncol = length(school_r), nrow = length(days), byrow = FALSE) > 
                         matrix(school_r, ncol = length(school_r), nrow = length(days), byrow = TRUE));

y_school <- data.frame(t = days, school = 1 - (n_closures - n_reopenings))
googmo <- merge(googmo, y_school, by = "t")
googmo[, date := ymd("2020-01-01") + t]

schedule <- make_schedule(googmo)

rm(y_school, googmo)

