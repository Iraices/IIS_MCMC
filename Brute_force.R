#### ### Brute force sampling in Stan
#### Meta analysis
###########################################################################

library(rstan)
library(readxl)
library(dplyr)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Load data
source('Load_biomanipulation_data.R')
## Load sets of priors
source('generate_sets_of_priors.R')

##################################################################################################
## Create functions to run the random effect linear model with a set of priors
## and find bounds on the expected effect
{
stan_meta_analysis = stan_model(file = "Meta_analysis_RE.stan")

brute_force_means = function(grid, y, sigma, n_iter = 10000){
  posterior_mu_means = rep(0, dim(grid)[1])
  tau_mu_means =  rep(0, dim(grid)[1])
  tau_study_means =  rep(0, dim(grid)[1])
  k_means =  rep(0, dim(grid)[1])
  log_k_means =  rep(0, dim(grid)[1])
  
   for(i in 1:dim(grid)[1]){
    model_meta_analysis = sampling(stan_meta_analysis, data = list(N=length(y), y = y, sigma = sigma, mu0 = grid$mu_interval[i], tau0 = grid$tau_interval[i]),
                                                           iter = n_iter, warmup = n_iter/2, chains = 2, seed = 5,
                                                            control = list(adapt_delta = 0.9, max_treedepth = 15))
    
    summary_model = summary(model_meta_analysis)
    
    posterior_mu_means[i] = summary_model$summary['mu','mean']
    tau_mu_means[i] = summary_model$summary['tau_mu','mean']
    tau_study_means[i] = summary_model$summary['tau_study','mean']
    k_means[i] = summary_model$summary['k','mean']
    log_k_means[i] = summary_model$summary['log_k','mean']
      
    print(i)
    }
  return(all_means = cbind(grid,  posterior_mu_means, tau_mu_means, tau_study_means, k_means, log_k_means))
}

## Minimun and maximum mean mu values
cal_min_max_means = function(output_data){
  
  lower_index = which.min(output_data$posterior_mu_means)
  
  upper_index = which.max(output_data$posterior_mu_means) 
  
  return(list(lower_mean = output_data$posterior_mu_means[lower_index], mu0_lower_mean = output_data$mu_interval[lower_index], tau0_lower_mean = output_data$tau_interval[lower_index], 
              upper_mean = output_data$posterior_mu_means[upper_index], mu0_upper_mean = output_data$mu_interval[upper_index], tau0_upper_mean = output_data$tau_interval[upper_index])
  )
}
}
###################################################################################################
#################################################################################################################
## Runs brute force for all available data (N = 75) 
{
initial_time_min = proc.time()

out_brute_force = brute_force_means(grid = grid, y = y, sigma = sigma, n_iter = 20000)

save(out_brute_force, file = "out_brute_force.Rdata")
load("out_brute_force.Rdata")

output = cal_min_max_means(out_brute_force)
 
final_time_min = proc.time() - initial_time_min 
}

## Study the results

output[,c("lower_mean","mu0_lower_mean","tau0_lower_mean")]

