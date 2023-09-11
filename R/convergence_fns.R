# Functions for convergence diagnostics

# detect outlier chains
outliers <- function(post_list, pops) {
  r = NULL;
  for (p in pops) {
    vars = names(post_list[[p]])[-1:-4]
    post = post_list[[p]]
    for (v in vars) {
      mu = post[, mean(get(v))]
      sd = post[, sd(get(v))]
      
      ind = post[, mean(get(v)), by = chain][, .(chain, div = abs(V1 - mean(V1)), indicator = abs(V1 - mu) > 1 * sd)]
      if (ind[, any(indicator == TRUE)]) {
        r = rbind(r, 
                  data.table(pop = p, variable = v, chains = ind[indicator == TRUE, chain])
        )
      }
    }
  }
  return (r)
}

# Function to remove outliers and burn-in

strip_burn <- function(post_list, pops, out, burn_in) {
  for (p in pops)
  {
    post_list[[p]][, trial := trial - burn_in]
    if (is.null(out)) {
      outchains = numeric()
    } else {
      outchains = out[pop == p, unique(chains)]
    }
    post_list[[p]] = post_list[[p]][trial >= 0 & !chain %in% outchains]
  }
  post_list
}

# Function for gelman-rubin diagnostic, Rhat

gelman_rubin <- function(post, varname) {
  # number of chains
  M = post[, uniqueN(chain)]
  
  # number of iterations
  N = post[, max(trial + 1), by = chain][, min(V1)]
  
  # posterior mean for each chain
  thm = post[, mean(get(varname)), by = chain][, V1]
  
  # posterior variance for each chain
  s2m = post[, var(get(varname)), by = chain][, V1]
  
  # grand posterior mean
  th = mean(thm)
  
  # how the individual means vary around the grand mean
  B_over_n = var(thm)
  
  # averaged variance of the chains
  W = mean(s2m)
  
  Vh = (N - 1) * W / N + (M + 1) * B_over_n / M;
  
  sqrt(Vh/W)
}

# Calculate all Rhats for a model fit
gr_all <- function(post_list, pops) {
  r = NULL
  for (p in pops) {
    vars = names(post_list[[p]])[-1:-4]
    for (v in vars) {
      g = gelman_rubin(post_list[[p]], v);
      r = rbind(r, data.table(pop = p, variable = v, Rhat = g))
    }
  }
  r
}

# Caterpillar plot for posteriors
caterpillar_plot <- function(post_list, pops)
{
  plots = list()
  for (p in pops) {
    pm = melt(post_list[[p]], id.vars = c("trial", "chain"))
    plots[[length(plots) + 1]] = ggplot(pm) +
      geom_line(aes(x = trial, y = value, group = chain)) +
      facet_wrap(~variable, scales = "free_y", ncol = 2) +
      theme_cowplot(font_size = 8) + labs(y = NULL, title = "Dominican Republic")
  }
  plot_grid(plotlist = plots, nrow = 1)
}

# Wrapper function for convergence diagnostics

convergence <- function(post_list, pop) {
  post <- post_list[[1]]
  cat("Detecting outlier chains...\n");
  out <- outliers(post, pop)
 
   cat("Calculating convergence diagnostics...\n");
  Rhat <- cbind(gr_all(post_list, pop))
  
  warn <- sum(Rhat$Rhat > 1.1, na.rm = TRUE)
  if(warn > 0) {print(paste0("Warning: ", warn, " Rhat > 1.1"))}
  cat("Building trace plot...\n");
  cp <- caterpillar_plot(post_list, pop)
  cat("Saving trace plot...\n");
  return(list(outliers = out, Rhat = Rhat, caterpillar = cp))
  
}

