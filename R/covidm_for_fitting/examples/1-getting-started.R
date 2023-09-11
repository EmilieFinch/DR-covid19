# 1-getting-started.R
# for covidm version 2

# covidm options
cm_path = "~/Projects/newcovid-master/covidm_for_fitting/"; ### CHANGE THIS to reflect the path to covidm.
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "R/covidm.R"))

# build parameters for all of UK, down to the regional level (level 2).
params = cm_parameters_SEI3R(cm_uk_locations("UK", 2), deterministic = T);

nr_pops <- 12
for(i in 1:nr_pops){
  params$pop[[i]]$dE2 <- params$pop[[i]]$dE
  params$pop[[i]]$dIp2 <- params$pop[[i]]$dIp
  params$pop[[i]]$dIa2 <- params$pop[[i]]$dIa
  params$pop[[i]]$dIs2 <- params$pop[[i]]$dIs
  params$pop[[i]]$seed_times2 <- params$pop[[i]]$seed_times2}


# # alternatively: build parameters for another country.
# params = cm_parameters_SEI3R("Italy");

# run the model
run = cm_simulate(params, 1)

# show results
ggplot(run$dynamics[compartment == "S"]) +
    geom_line(aes(t, value, colour = group, group = group)) +
    facet_wrap(~population)

