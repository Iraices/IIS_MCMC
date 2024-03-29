## load ranges for sets of priors
load('setofpriorstouse.Rdata')

low_mu0 <- mu0_set_prior[1]
high_mu0 <- mu0_set_prior[2]
low_tau0 <- tau0_set_prior[1] 
high_tau0 <- tau0_set_prior[2]

n_iter_stan <- 20000
mu_interval  <- seq(from = low_mu0, to = high_mu0, by = 0.5)
tau_interval <- seq(from = low_tau0, to = high_tau0, by = 0.25)
grid <- expand.grid(mu_interval = mu_interval, tau_interval = tau_interval)
