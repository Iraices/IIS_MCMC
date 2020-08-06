### Load model
stan_meta_analysis = stan_model(file = "Meta_analysis_RE.stan")

## Functions
  
  get_log_posterior_MCMC_sample = function(mu0, log_tau0, y, sigma, n_iter_stan, target_ESS){
    tau0  = exp(log_tau0)  
    warmup = round(0.5 * n_iter_stan)
    
    repeat{
      model_meta_analysis = sampling(stan_meta_analysis, 
                                     data = list(N = length(y), y = y, sigma = sigma, mu0 = mu0, tau0 = tau0),
                                     iter = n_iter_stan, warmup = warmup, chains = 2, seed = 5, 
                                     control = list(adapt_delta = 0.9, max_treedepth = 15))
      
      postsam = rstan::extract(model_meta_analysis, permuted = FALSE) 
      
      mu = c(postsam[,1,'mu'],postsam[,2,'mu'])
      tau_mu = c(postsam[,1,'tau_mu'], postsam[,2,'tau_mu'])
      tau_study = c(postsam[,1,'tau_study'], postsam[,2,'tau_study'])
      k = c(postsam[,1,'k'], postsam[,2,'k'])
      log_k = c(postsam[,1,'log_k'], postsam[,2,'log_k'])
      log_posterior = c(postsam[,1,'log_posterior'], postsam[,2,'log_posterior'])
      
      # ESS_MCMC is the effective sample size for MCMC samples
      cor = acf(mu, plot = FALSE)
      index = which(cor$acf < 0)[1]
      if(index==2){
        ESS_MCMC_mu = length(mu)  
      }
      else{
        ESS_MCMC_mu = length(mu) / (1 + 2 * sum(cor$acf[2:(index - 1)]))
      }
      
      ## create init function in case of a low sample size
      if(ESS_MCMC_mu < 1.2 * target_ESS){ # target_ESS + 20% target_ESS
        n_iter_stan = n_iter_stan + 10000
        warmup = round(0.5 * n_iter_stan)
      }
      else{
        break
      }
    }
    return(list(log_posterior = log_posterior,
                post = cbind(mu = mu, tau_mu = tau_mu, k = k, log_k = log_k, tau_study = tau_study),  
                ESS_MCMC_mu =  ESS_MCMC_mu,  hyp_par = c(mu0 = mu0, log_tau0 = log_tau0)))
  }
  
  get_log_posterior_target_sample = function(mu0, log_tau0, log_posterior_MCMC_sample, y, sigma){
    
    mu = log_posterior_MCMC_sample$post[,'mu']
    tau_mu = log_posterior_MCMC_sample$post[,'tau_mu']
    tau_study = log_posterior_MCMC_sample$post[,'tau_study']
    k = log_posterior_MCMC_sample$post[,'k']
    log_k = log_posterior_MCMC_sample$post[,'log_k']
    
    tau0 = exp(log_tau0)
    
    N = length(y)
    S = length(mu)
    
    log_lik_data = matrix(0, S, N)
    
    for (i in 1:N){
      log_lik_data[,i] = dnorm(y[i], mean =  mu, sd = sqrt(sigma[i] * sigma[i] + tau_study * tau_study), log = TRUE)
    }
    log_lik = rowSums(log_lik_data)
    log_posterior_target_sample = log_lik + dnorm(mu, mean = mu0, sd = tau_mu, log = TRUE) + dnorm(tau_mu, mean = 0, sd = tau0,log = TRUE) + dnorm(log_k, mean = 0, sd = 1.5, log = TRUE)
    return(list(log_likelihood = log_lik, log_posterior = log_posterior_target_sample))
  }
  
  calc_weight = function(log_posterior_MCMC_sample, log_posterior_target_sample, mu_posterior_MCMC_sample){
    
    number_samples = length(mu_posterior_MCMC_sample)
    
    log_weight = log_posterior_target_sample - log_posterior_MCMC_sample
    # weight estimate
    weight = exp(log_weight)
    # importance sampling expectation
    mu_x_weight = mu_posterior_MCMC_sample * weight
    weight_sum = sum(weight)
    
    IS_Est_mu = sum(mu_x_weight) / weight_sum
    
    # ESS_MCMC is the effective sample size for MCMC samples
    g = mu_x_weight - (IS_Est_mu * weight)
    cor = acf(g, plot = FALSE)
    index = which(cor$acf < 0.05)[1]
    if(index==2){
    ESS_MCMC = number_samples  
    }
    else{
    ESS_MCMC = number_samples / (1 + 2 * sum(cor$acf[2:(index - 1)]))
    }
    #ESS_IS is the effective sample size for standard IS 
    ESS_IS = (weight_sum^2) / (sum(weight^2))
    
    # Total effective sample size ESS
    ESS = (ESS_MCMC / number_samples) * ESS_IS
    
    return(list(weight = weight, IS_Est_mu = IS_Est_mu, g = g, 
                ESS_IS = ESS_IS, ESS_MCMC = ESS_MCMC, ESS = ESS))
  }
  
  calc_EV_with_penalty = function(parameters, lower_bound_mu0, upper_bound_mu0, lower_bound_tau0, upper_bound_tau0, 
                                  log_posterior_MCMC_sample, y, sigma, y_int, region){
    
    penalty_term = 1000*max(abs(y))
    mu0 = parameters[1]
    log_tau0 = parameters[2]
    
    log_post_target_sample = get_log_posterior_target_sample(mu0 = mu0, log_tau0 = log_tau0, 
                              log_posterior_MCMC_sample = log_posterior_MCMC_sample, y = y, sigma = sigma)$log_posterior
    
    w = calc_weight(log_posterior_MCMC_sample = log_posterior_MCMC_sample$log_posterior, 
                    log_posterior_target_sample = log_post_target_sample, 
                    mu_posterior_MCMC_sample = log_posterior_MCMC_sample$post[,'mu'])
    
    ## Introduce a penalty if the choice of prior falls outside the set of priors
    ## The penalty term is set to be larger than the IS_Est_mu (e.g.ensure it is larger than the sample mean of data)
    ## The if condition checks that selected priors satisfice the  requirements of Prior Predictive Check
    if(mu0 >= lower_bound_mu0 & mu0 <=  upper_bound_mu0 & log_tau0 >=lower_bound_tau0 & log_tau0 <= upper_bound_tau0 ){
      
      eval = exp(log_tau0) - region(mu0)
      
      estimate_mean_with_penalty = w$IS_Est_mu + penalty_term * max(0, eval)
    }
    else{
      estimate_mean_with_penalty = 2*penalty_term 
    }
    return(estimate_mean_with_penalty)
  } 
  
  estimate_lower_bound = function(mu0, log_tau0, lower_bound_mu0, upper_bound_mu0, lower_bound_tau0, upper_bound_tau0, 
                                  y, sigma, y_int, n_iter_stan, region, target_ESS){
    
    number_of_iteration = 0
    
    par_init = c(mu0, log_tau0)
    
    save_status <- c()
    
    repeat{
      
      log_posterior_MCMC_sample = get_log_posterior_MCMC_sample(mu0 = mu0, log_tau0 = log_tau0, y = y, sigma = sigma, n_iter_stan = n_iter_stan, target_ESS = target_ESS)
      
      number_samples = length(log_posterior_MCMC_sample$post[,'mu'])
      
      ## optim function gives as outputs the priors that minimized the expected value and satisfice the requiriments 
      ## of the Prior Predictive Check 
      argmin = optim(par = par_init, fn = calc_EV_with_penalty, method = "SANN", 
                     #control = list(maxit = 20000, tmax = 20),
                     lower_bound_mu0 = lower_bound_mu0, upper_bound_mu0 = upper_bound_mu0, 
                     lower_bound_tau0 = lower_bound_tau0, upper_bound_tau0 = upper_bound_tau0,
                     log_posterior_MCMC_sample = log_posterior_MCMC_sample, y = y, sigma = sigma, y_int = y_int,
                     region = region)
      
      mu0 = argmin$par[1]
      log_tau0 = argmin$par[2]
      par_init = c(mu0, log_tau0)
      
      log_posterior_target_sample =  get_log_posterior_target_sample(mu0 =  mu0, log_tau0 = log_tau0, 
                                                                     log_posterior_MCMC_sample = log_posterior_MCMC_sample, 
                                                                     y = y, sigma = sigma)
      
      w = calc_weight(log_posterior_MCMC_sample = log_posterior_MCMC_sample$log_posterior, 
                      log_posterior_target_sample = log_posterior_target_sample$log_posterior, 
                      mu_posterior_MCMC_sample = log_posterior_MCMC_sample$post[,'mu'])
      
      number_of_iteration = number_of_iteration + 1
      
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
  