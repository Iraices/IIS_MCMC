#### Iterative Importance sampling under sets of priors using MCMC sampling
###########################################################################

{
## Load data
source('Load_biomanipulation_data.R')
  
data = data.frame(y = y, sigma = sigma)

## Generate data sets with smaller sample sizes

random_data = sample_n(data, size = 40, replace = FALSE)

random_data_precise = data[order(data$sigma),][1:40,]

save(random_data, file = "random_data.Rdata")
save(random_data_precise, file = "random_data_precise.Rdata")

load("random_data.Rdata")
load("random_data_precise.Rdata")

}
#######
## Load sets of priors
source('generate_sets_of_priors.R')

## Create functions
source('Create_Functions_For_IIS.R')

## IIS based on full data set with changed estimation errors
{

initial_time = proc.time()

output_orig = estimate_lower_bound (mu0 = 0, log_tau0 = log(5), 
                                     lower_bound_mu0 = low_mu0, upper_bound_mu0 = high_mu0, 
                                     lower_bound_tau0 = low_log_tau0, upper_bound_tau0 = high_log_tau0, 
                                     y = data$y, sigma = data$sigma, y_int = y_int,
                                     n_iter_stan = n_iter_stan, region = reg, target_ESS = 5000)

final_time = proc.time() - initial_time

save(output_orig, file = "output_orig.Rdata")
load("output_orig.Rdata")

###################################################

initial_random_data_time = proc.time()

output_random_data = estimate_lower_bound (mu0 = 0, log_tau0 = log(5), 
                                    lower_bound_mu0 = low_mu0, upper_bound_mu0 = high_mu0, 
                                    lower_bound_tau0 = low_log_tau0, upper_bound_tau0 = high_log_tau0, 
                                    y = random_data$y, sigma = random_data$sigma, y_int = y_int,
                                    n_iter_stan = n_iter_stan, region = reg, target_ESS = 5000)

final_random_data_time = proc.time() - initial_random_data_time

save(output_random_data, file = "output_random_data.Rdata")
load("output_random_data.Rdata")

##################################################

initial_random_data_precise_time = proc.time()

output_random_data_precise = estimate_lower_bound (mu0 = 0, log_tau0 = log(5), 
                                    lower_bound_mu0 = low_mu0, upper_bound_mu0 = high_mu0, 
                                    lower_bound_tau0 = low_log_tau0, upper_bound_tau0 = high_log_tau0, 
                                    y = random_data_precise$y, sigma = random_data_precise$sigma, y_int = y_int,
                                    n_iter_stan = n_iter_stan, region = reg, target_ESS = 5000)

final_random_data_precise_time = proc.time() - initial_random_data_precise_time

save(output_random_data_precise, file = "output_random_data_precise.Rdata")
load("output_random_data_precise.Rdata")

}
