#### Iterative Importance sampling under sets of priors using MCMC sampling
###########################################################################

{
## Load data
source('Load_biomanipulation_data.R')
  
data <- data.frame(y = y, sigma = sigma)

}
#######
## Load sets of priors
source('generate_sets_of_priors.R')

## Create functions
source('Functions_For_IIS.R')

## IIS based on different target_ESS and initial values 
{

output_orig_1 <- estimate_lower_bound (mu0 = -7, tau0 = 6, 
                                     lower_bound_mu0 = low_mu0, upper_bound_mu0 = high_mu0, 
                                     lower_bound_tau0 = low_tau0, upper_bound_tau0 = high_tau0, 
                                     y = data$y, sigma = data$sigma, y_int = y_int,
                                     n_iter_stan = n_iter_stan, region = reg, target_ESS = 5000)

save(output_orig_1, file = "output_orig_1.Rdata")
load("output_orig_1.Rdata")

###################################

output_orig_2 <- estimate_lower_bound (mu0 = -7, tau0 = 6, 
                                       lower_bound_mu0 = low_mu0, upper_bound_mu0 = high_mu0, 
                                       lower_bound_tau0 = low_tau0, upper_bound_tau0 = high_tau0, 
                                       y = data$y, sigma = data$sigma, y_int = y_int,
                                       n_iter_stan = n_iter_stan, region = reg, target_ESS = 10000)

save(output_orig_2, file = "output_orig_2.Rdata")
load("output_orig_2.Rdata")

###################################

output_orig_3 <- estimate_lower_bound (mu0 = 10, tau0 = 10, 
                                       lower_bound_mu0 = low_mu0, upper_bound_mu0 = high_mu0, 
                                       lower_bound_tau0 = low_tau0, upper_bound_tau0 = high_tau0, 
                                       y = data$y, sigma = data$sigma, y_int = y_int,
                                       n_iter_stan = n_iter_stan, region = reg, target_ESS = 10000)

save(output_orig_3, file = "output_orig_3.Rdata")
load("output_orig_3.Rdata")

}

