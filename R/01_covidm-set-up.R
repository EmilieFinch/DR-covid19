# -------------------------------------------------------------------------------------------------------
# Set up covidM for the Dominican Republic
# Author: Emilie Finch, adapted from code by Nick Davies in https://github.com/nicholasdavies/newcovid
# -------------------------------------------------------------------------------------------------------

# Source fn files  --------------------------------------------------------------------------------------

source(here("R", "cpp_fns.R"))
source(here("R", "run-covidM_fns.R"))
source(here("R", "generate-fit_fns.R"))
source(here("R", "convergence_fns.R"))
source(here("R", "counterfactual_fns.R"))

# Load covidM

cm_path <- here("R", "covidm_for_fitting")
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(here(cm_path, "R/covidm.R"))


# Set number of iterations for burn-in and MCMC --------------------------------------------------------

REP_START = 1
REP_END = 1
BURN_IN = 3000
ITER = 5000
BURN_IN_FINAL = 3000
ITER_FINAL = 5000 
N_THREADS = 10; # Number of threads for fitting, set to number of cores 


# Set options for covidM run ---------------------------------------------------------------------------

# These are set to the options used for the main manuscript results
# With other options for sensitivity analyses in comments

date_fitting <- "2022-03-01" # Date until which model simulations are produced
opt_mobility <- "comix-base" # Mobility assumption, e.g. here a scaling parameter is fitted between comix-adjusted mobility 
                             # and baseline mobility, can be one of: comix-base, comix-only
opt_wane <- "central" # Waning assumption, can be one of: central, high, low or none
seas_yn <- "seasno" # Include seasonality or not, can be: seasyes or seasno
seas_amp <- 0.2 # If including seasonality, amplitude. Not used if seas_yn set to "seasno"
extra_voc_takeoff <- "central" # Timing for Mu transmission increase, can be: central, broad or late
opt_v2 = TRUE # Introduce second variant: TRUE/FALSE
opt_relu = TRUE # Fit parameter to modify transmissibility of second variant: TRUE/FALSE

# Build parameters for Dominican Republic --------------------------------------------------------------

which_pops = 1 # Number of regions to fit

params = cm_parameters_SEI3R("Dominican Republic", deterministic = T, 
                             date_start = "2020-01-01", 
                             date_end = date_fitting,
                             dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
                             dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
                             dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
                             dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p)
params = cm_split_matrices_ex_in(params, 15)

params$pop[[1]]$rho = rep(1, 16)
params$schedule = schedule

# Load age-specific susceptibility and probability of developing symptoms-------------------------------

covid_rates = readRDS(here("data", "age-susceptibility-symptomatic.rds"))
covu = covid_rates$covu;
covy = covid_rates$covy;

for (i in seq_along(params$pop)) {
  params$pop[[i]]$u = covu / mean(covu);
  params$pop[[i]]$u2 = covu / mean(covu);
  params$pop[[i]]$y = covy;
  params$pop[[i]]$y2 = covy;
}

# Vaccine efficacy parameters --------------------------------------------------------------------------

# Note ei_v and ed_v1 refer to strain 1 and ei2_v and ed_vi2 refer to strain 2.

params$pop[[1]]$ei_v = rep(0.67,16) # Efficacy against infection for WT / Alpha and other intermediate VOI / VOC
ed_oi = 0.67 # Overall efficacy against disease for WT / Alpha and other intermediate VOI/VOC
params$pop[[1]]$ei2_v = rep(0.39 ,16) # Efficacy against infection for Delta
ed_oi2 = 0.39 # Overall efficacy against disease for Delta
params$pop[[1]]$ed_vi = calc_ve_d(params$pop[[1]]$ei_v, ed_oi) # Efficacy against disease given infection, that is clinical fraction among breakthrough for WT / Alpha and other intermediate VOI / VOC
params$pop[[1]]$ed_vi2 = calc_ve_d(params$pop[[1]]$ei2_v, ed_oi2) # Efficacy against disease given infection, that is clinical fraction among breakthrough for Delta

# Waning parameters -------------------------------------------------------------------------------------

if(opt_wane == "central"){
  # using assumption of 15% of individuals with immunity having reinfections within a year
  params$pop[[1]]$wn = rep(log(0.85)/-365,16)
  params$pop[[1]]$wn2 = rep(log(0.85)/-365, 16)
  params$pop[[1]]$wv = rep(log(0.6)/-182.5, 16) # Assuming 40% of individuals will have reinfections within 6 months. Based on Cerqueira-Silva et al 2022.
  
}

if(opt_wane == "high"){
  params$pop[[1]]$wn = rep(log(0.85)/-182.5,16)
  params$pop[[1]]$wn2 = rep(log(0.85)/-182.5, 16)
  params$pop[[1]]$wv = rep(log(0.6)/-91.25, 16) # Assuming 40% of individuals will have reinfections within 3 months. Based on Cerqueira-Silva et al 2022.
  
}

if(opt_wane == "low"){
  params$pop[[1]]$wn = rep(log(0.85)/-730,16)
  params$pop[[1]]$wn2 = rep(log(0.85)/-730, 16)
  params$pop[[1]]$wv = rep(log(0.84)/-120, 16) # Assuming 16% of individuals will have reinfections within 4 months. Based on Cerqueira-Silva et al 2022 Supp Table 4.
  
}

if(opt_wane == "none"){
  params$pop[[1]]$wn = rep(0,16)
  params$pop[[1]]$wn2 = rep(0, 16)
  params$pop[[1]]$wv = rep(0, 16)
  
}

# Transmission increase for extra voc --------------------------------------------------------------------------

if(extra_voc_takeoff == "broad"){
  voct <- seq(404,495, by = 1) # 8th Feb to 10th May
} else if(extra_voc_takeoff == "central"){
  voct <- seq(404,460, by = 1) # 8th Feb to 5th April
} else if(extra_voc_takeoff == "late"){
  voct <- seq(432,460, by = 1) # 8th March to 5th April
} else {}

voci <- seq(0,1, by = 1/(length(voct) -1))
