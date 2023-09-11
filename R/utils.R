
# Function to create schedule for covidm input from Google data and school closures
make_schedule <- function(y) {

  schedule = list()
  schedule[[1]]  = list(parameter = "contact", pops = 0,  mode = "assign", values = list(), times = numeric(0))
  
  for (r in 1:nrow(y)){
    schedule[[1]]$values = c(schedule[[1]]$values,
                             list(y[r, c(residential, workplace, grocpharm, retrec, transit, school)]));
    schedule[[1]]$times = c(schedule[[1]]$times, y[r, t])
  }
  
  return (schedule)
}

# Function to read in immune states output

read_dynamics <- function(){
  files <- list.files(here("output", "model-dynamics"))
  dyn_out <- data.frame()
  for (f in 1:length(files)){
    temp <- qread(here("output", "model-dynamics", files[f]))
    dyn_out <- rbind(dyn_out, temp)
  }
  return(dyn_out)
}

# Distribution functions

skewness <- function(x) {
  mean((x - mean(x))^3) / (mean((x - mean(x))^2)^1.5)
}

skew_normal_objective = function(x, target) {
  central = target[1]
  lower = target[2]
  upper = target[3]
  
  xi = x[1]
  omega = x[2]
  alpha = x[3]
  
  delta = alpha / sqrt(1 + alpha^2)
  
  mean = xi + omega * delta * sqrt(2 / pi)
  lq = sn::psn(lower, xi, omega, alpha)
  uq = sn::psn(upper, xi, omega, alpha)
  
  return ((central - mean) ^ 2 + (lq - 0.025) ^ 2 + (uq - 0.975) ^ 2)
}

skew_normal_solve <- function(central, lower, upper, adj_react2)
{
  d = NULL;
  set.seed(12345)
  for (i in seq_along(central))
  {
    if (adj_react2[i]) {
      # Get point estimate of q, crude seroprevalence (non-adjusted)
      p = central[i]
      sensitivity = 0.844
      specificity = 0.986
      q = p * (sensitivity) + (p - 1) * (specificity - 1)
      
      # Model p with uncertainty over sensitivity and specificity
      sensitivity = rbeta(100000, shape1 = 31.65, shape2 = 6.3)
      specificity = rbeta(100000, 4*98.6, 4*1.4)[sensitivity > q]
      sensitivity = sensitivity[sensitivity > q]
      p = (q + specificity - 1) / (sensitivity + specificity - 1);
      
      # Get skew-normal distribution that matches these parameters
      par = cp2dp(c(mean(p), sd(p), min(0.98, skewness(p))), "SN");
      d = rbind(d,
                data.table(xi = par[1], omega = par[2], alpha = par[3]));
    } else {
      solution = optim(par = c(central[i], (upper[i] - lower[i]) / 4, 0), fn = skew_normal_objective, 
                       target = c(central[i], lower[i], upper[i]), 
                       method = "Nelder-Mead", control = list(maxit = 10000, abstol = 0, reltol = 0));
      d = rbind(d, 
                data.table(xi = solution$par[1], omega = solution$par[2], alpha = solution$par[3])
      );
    }
  }
  
  return (d)
}

# Function for normally distributed delay

delay_normal <- function(mu, sd, t_max, t_step)
{
  t_points = seq(0, t_max, by = t_step);
  heights = pnorm(t_points + t_step/2, mu, sd) -
    pnorm(pmax(0, t_points - t_step/2), mu, sd);
  return (data.table(t = t_points, p = heights / sum(heights)))
}

# Function to calculate efficacy against disease given infection from overall efficacy against disease and efficacy against infection

calc_ve_d = function(ve_i, ve_d_o){
  ve_d = (ve_d_o - ve_i)/ (1 - ve_i)
  
  return(ve_d)
}

pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100


# Helper function to generate changes in contacts over time

generate_contacts <- function(schedule, posteriorsI, mobility_type){
  
  # Create table of mobility plot 
  # Curves from Co-Mix analysis
  
  work_curve <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361) 
  other_curve <- c(0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701)
  
  wplc <- c()
  groc <- c()
  rtrc <- c()
  trns <- c()
  scho <- c()
  othx <- c()
  bwplc <- c()
  bothx <- c()
  t <- c()
  home  <- c()
  
  for (i in 1:length(schedule[[1]]$values)) {
    home[i] = schedule[[1]]$values[[i]][[1]]
    wplc[i] = schedule[[1]]$values[[i]][2]
    groc[i] = schedule[[1]]$values[[i]][3]
    rtrc[i] = schedule[[1]]$values[[i]][4]
    trns[i] = schedule[[1]]$values[[i]][5]
    scho[i] = schedule[[1]]$values[[i]][6]
    othx[i] = min((as.numeric(rtrc[i]) * 0.345 + as.numeric(trns[i]) * 0.445 + as.numeric(groc[i]) * 0.210), 1.24)
    bwplc[i] = schedule[[1]]$values[[1]][[2]] # Note this is 1
    bothx[i] = schedule[[1]]$values[[1]][[3]] # Note this is 1
    t[i] = schedule[[1]]$times[i]
  }
  
  # Create vectors for baseline google mobility
  
  cppFunction(
    'auto interp(double y, NumericVector curve) {
        unsigned int i = (unsigned int)(y * 100);
        double f = y * 100 - i;
        if (y < 0) {
        return curve[0];
        } else if (i >= curve.size() - 1){
         return curve[126];
        } else {
        return (1 - f) * curve[i] + f * curve[i + 1];
        }
    }')
  
  # Modify based on posterior values for contact adjustment 
  contact_adj_a <- density(posteriorsI[[1]]$contact_adj_a)$x[which(density(posteriorsI[[1]]$contact_adj_a)$y == max(density(posteriorsI[[1]]$contact_adj_a)$y))]
  contact_adj_b <- density(posteriorsI[[1]]$contact_adj_b)$x[which(density(posteriorsI[[1]]$contact_adj_b)$y == max(density(posteriorsI[[1]]$contact_adj_b)$y))]
  
  new_wplc <- c()
  new_othx <- c()

  if(mobility_type == "comix-base") {
    for(i in 1:365){
      new_wplc[i] <- min((contact_adj_a*interp(wplc[i], work_curve) + ((1-contact_adj_a)*bwplc[i])), 1.0)
      new_othx[i] <- min((contact_adj_a*interp(othx[i], other_curve) + ((1-contact_adj_a)*bothx[i])), 1.0)
    }  
    for(i in 366:length(wplc)){
      new_wplc[i] <- min((contact_adj_b*interp(wplc[i], work_curve) + ((1-contact_adj_b)*bwplc[i])), 1.0)
      new_othx[i] <- min((contact_adj_b*interp(othx[i], other_curve) + ((1-contact_adj_b)*bothx[i])), 1.0)
    } 
  } else if(mobility_type == "comix-only" ){
    
    for(i in 1:length(wplc)){
      new_wplc[i] <- min(interp(wplc[i], work_curve), 1.0)
      new_othx[i] <- min(interp(othx[i], other_curve), 1.0)
    }
  }
  comix_wplc <- c()
  comix_othx <- c()
  for(i in 1:length(wplc)) {
    comix_wplc[i] <- interp(wplc[i], work_curve)
    comix_othx[i] <- interp(othx[i], other_curve)
  }
  
  fin_wplc <- new_wplc
  fin_othx <- new_othx
  fin_home <- home
  fin_scho <- scho
  
  return(list(contact_times = t, fin_wplc = fin_wplc, fin_othx = fin_othx, fin_home = fin_home, fin_scho = fin_scho, comix_wplc = comix_wplc, comix_othx = comix_othx, wplc = wplc, othx = othx))
}

