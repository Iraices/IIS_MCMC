### Load model
stan_meta_analysis = stan_model(file = "Meta_analysis_RE.stan")

## Functions

get_log_posterior_MCMC_sample <- function(mu0_MCMC_sample, tau0_MCMC_sample, y, sigma, n_iter_stan, target_ESS){
  warmup <- round(0.5 * n_iter_stan)
  
  repeat{
    model_meta_analysis <- sampling(stan_meta_analysis, 
                                    data = list(N = length(y), y = y, sigma = sigma, mu0 = mu0_MCMC_sample, tau0 = tau0_MCMC_sample),
                                    iter = n_iter_stan, warmup = warmup, chains = 2, seed = 5, 
                                    control = list(adapt_delta = 0.9, max_treedepth = 15))
    
    postsam <- rstan::extract(model_meta_analysis, permuted = FALSE) 
    
    mu <- c(postsam[,1,'mu'],postsam[,2,'mu'])
    tau_mu <- c(postsam[,1,'tau_mu'], postsam[,2,'tau_mu'])
    tau_study <- c(postsam[,1,'tau_study'], postsam[,2,'tau_study'])
    k <- c(postsam[,1,'k'], postsam[,2,'k'])
    posterior <- exp(c(postsam[,1,'log_posterior'], postsam[,2,'log_posterior']))
    
    ESS_MCMC_mu <- summary(model_meta_analysis)$summary['mu','n_eff']
    
    ## in case of a low sample size
    if(ESS_MCMC_mu < 1.2 * target_ESS){ # target_ESS + 20% target_ESS
      n_iter_stan <- n_iter_stan + 20000
      warmup <- round(0.5 * n_iter_stan)
    }
    else{
      break
    }
  }
  return(list(posterior = posterior,
              post = cbind(mu = mu, tau_mu = tau_mu, k = k, tau_study = tau_study),  
              ESS_MCMC_mu =  ESS_MCMC_mu,  hyp_par = c(mu0_MCMC_sample = mu0_MCMC_sample, tau0_MCMC_sample = tau0_MCMC_sample)))
}

get_log_posterior_target_sample <- function(mu0_target_sample, tau0_target_sample, log_posterior_MCMC_sample, y, sigma){
  
  mu <- log_posterior_MCMC_sample$post[,'mu']
  tau_mu <- log_posterior_MCMC_sample$post[,'tau_mu']
  tau_study <- log_posterior_MCMC_sample$post[,'tau_study']
  k <- log_posterior_MCMC_sample$post[,'k']
  
  N <- length(y)
  S <- length(mu)
  
  log_lik_data <- matrix(0, S, N)
  
  
  for (i in 1:N){
    log_lik_data[,i] <- dnorm(y[i], mean =  mu, sd = sqrt(sigma[i] * sigma[i] + tau_study * tau_study), log = TRUE)
  }
  log_lik <- rowSums(log_lik_data)
  lik <- exp(log_lik)
  
  posterior_target_sample <- lik * dnorm(mu, mean = mu0_target_sample, sd = tau_mu) * dunif(tau_mu, min = 1, max = tau0_target_sample) * dunif(k, min = 1, max = 5)
  
  return(list(likelihood = lik, posterior = posterior_target_sample, 
              hyp_par = c(mu0_target_sample = mu0_target_sample, tau0_target_sample = tau0_target_sample)))
}

calc_weight <- function(posterior_MCMC_sample, posterior_target_sample, mu_posterior_MCMC_sample){
  
  number_samples <- length(mu_posterior_MCMC_sample)
  
  # weight estimate
  weight <- posterior_target_sample /  posterior_MCMC_sample
  
  # importance sampling expectation
  mu_x_weight <- mu_posterior_MCMC_sample * weight
  weight_sum <- sum(weight)
  
  IS_Est_mu <- sum(mu_x_weight) / weight_sum
  
  # ESS_MCMC is the effective sample size for MCMC samples
  g <- mu_x_weight - (IS_Est_mu * weight)
  cor <- acf(g, plot = FALSE)
  index <- which(cor$acf < 0.01)[1]
  if(is.na(index)){
    ESS_MCMC <- 1 
  }
  else if(index == 2){
    ESS_MCMC <- number_samples  
  }
  else{
    ESS_MCMC <- number_samples / (1 + 2 * sum(cor$acf[2:(index - 1)]))
  }
  #ESS_IS is the effective sample size for standard IS 
  ESS_IS <- (weight_sum^2) / (sum(weight^2))
  
  # Total effective sample size ESS
  ESS <- (ESS_MCMC / number_samples) * ESS_IS
  
  return(list(weight = weight, IS_Est_mu = IS_Est_mu, g = g, 
              ESS_IS = ESS_IS, ESS_MCMC = ESS_MCMC, ESS = ESS))
}

calc_EV_with_penalty <- function(parameters, lower_bound_mu0, upper_bound_mu0, lower_bound_tau0, upper_bound_tau0, 
                                 log_posterior_MCMC_sample, y, sigma, y_int, region){
  
  mu0 <- parameters[1]
  tau0 <- parameters[2]
  
  ## Introduce a penalty if the choice of prior falls outside the set of priors
  ## The penalty term is set to be larger than the IS_Est_mu (e.g.ensure it is larger than the sample mean of data)
  ## The if condition checks that selected priors satisfy the  requirements of Prior Predictive Check
  
  if(mu0 >= lower_bound_mu0 & mu0 <=  upper_bound_mu0 & tau0 >= lower_bound_tau0 & tau0 <= upper_bound_tau0){
    
    post_target_sample <- get_log_posterior_target_sample(mu0_target_sample = mu0, tau0_target_sample = tau0, 
                                                          log_posterior_MCMC_sample = log_posterior_MCMC_sample, y = y, sigma = sigma)$posterior
    
    w = calc_weight(posterior_MCMC_sample = log_posterior_MCMC_sample$posterior, 
                    posterior_target_sample = post_target_sample, 
                    mu_posterior_MCMC_sample = log_posterior_MCMC_sample$post[,'mu'])
    
    eval = tau0 - region(mu0)
    
    return(w$IS_Est_mu + penalty_term * max(0, eval))
  }
  return(1000*max(abs(y)))
}

estimate_lower_bound <- function(mu0, tau0, lower_bound_mu0, upper_bound_mu0, lower_bound_tau0, upper_bound_tau0, 
                                 y, sigma, y_int, region, n_iter_stan, target_ESS){
  
  number_of_iteration <- 0
  
  par_init <- c(mu0, tau0)
  
  save_status <- c()
  
  repeat{
    
    tau0_MCMC_sample <- tau0
    log_posterior_MCMC_sample <- get_log_posterior_MCMC_sample(mu0_MCMC_sample = mu0, tau0_MCMC_sample = tau0_MCMC_sample, 
                                                               y = y, sigma = sigma, n_iter_stan = n_iter_stan, target_ESS = target_ESS)
    
    number_samples <- length(log_posterior_MCMC_sample$post[,'mu'])
    
    argmin <- optim(par = par_init, fn = calc_EV_with_penalty, method = 'SANN',
                    lower_bound_mu0 = lower_bound_mu0, upper_bound_mu0 = upper_bound_mu0, 
                    lower_bound_tau0 = lower_bound_tau0, upper_bound_tau0 = upper_bound_tau0,
                    log_posterior_MCMC_sample = log_posterior_MCMC_sample, y = y, sigma = sigma, y_int = y_int,
                    region = region)
    
    mu0 <- argmin$par[1]
    tau0 <- argmin$par[2]
    par_init <- c(mu0, tau0)
    
    log_posterior_target_sample <- get_log_posterior_target_sample(mu0_target_sample = mu0, tau0_target_sample = tau0, 
                                                                   log_posterior_MCMC_sample = log_posterior_MCMC_sample, 
                                                                   y = y, sigma = sigma)
    
    w <- calc_weight(posterior_MCMC_sample = log_posterior_MCMC_sample$posterior, 
                     posterior_target_sample = log_posterior_target_sample$posterior, 
                     mu_posterior_MCMC_sample = log_posterior_MCMC_sample$post[,'mu'])
    
    number_of_iteration <- number_of_iteration + 1
    
    save_status <- rbind(save_status,
                         c(iteration_number = number_of_iteration, number_samples = number_samples,  mu_est = w$IS_Est_mu, ESS_IS = w$ESS_IS, ESS_MCMC = w$ESS_MCMC, ESS_MCMC_mu = log_posterior_MCMC_sample$ESS_MCMC_mu, ESS = w$ESS, par = par_init))
    
    print(list(iteration_number = number_of_iteration, number_samples = number_samples, mu_est = w$IS_Est_mu, ESS_IS = w$ESS_IS, ESS_MCMC = w$ESS_MCMC, ESS_MCMC_mu = log_posterior_MCMC_sample$ESS_MCMC_mu, ESS = w$ESS, par = par_init))
    
    if(w$ESS > target_ESS  || number_of_iteration == 10000){ 
      break  
    }
  }
  print(save_status)
  return(save_status)
}

