### Test function

library(rstan)
library(readxl)
library(dplyr)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Load data
source('Load_biomanipulation_data.R')

## Load sets of priors
source('generate_sets_of_priors.R')

## Create functions for IIS
source('Create_Functions_For_IIS.R')

########################################################
## Function

all_target_sample_weight_given_a_MCMC_sample = function(low_mu0, high_mu0, low_tau0, high_tau0, L, y, sigma, 
                                                        n_iter_stan, posterior_MCMC){
  
  
  mu_interval  = seq(from = low_mu0, to = high_mu0, length.out = L[1])
  tau_interval = seq(from = low_tau0, to = high_tau0, length.out = L[2])
  
  IS_Est_mu = matrix(0, nrow = L[2], ncol = L[1])
  ESS_IS = matrix(0, nrow = L[2], ncol = L[1])
  ESS_MCMC = matrix(0, nrow = L[2], ncol = L[1])
  ESS = matrix(0, nrow = L[2], ncol = L[1])
  rownames(IS_Est_mu) =  exp(tau_interval)
  colnames(IS_Est_mu) =  mu_interval
  rownames(ESS_IS) =  exp(tau_interval)
  colnames(ESS_IS) =  mu_interval
  rownames(ESS_MCMC) =  exp(tau_interval)
  colnames(ESS_MCMC) =  mu_interval
  rownames(ESS) =  exp(tau_interval)
  colnames(ESS) =  mu_interval
    
    for (i in 1:L[2]){ # tau0
      for (j in 1:L[1]){ # mu0
        log_posterior_target_sample = get_log_posterior_target_sample(mu0 = mu_interval[j], tau0 = tau_interval[i],
                                                                 log_posterior_MCMC_sample = posterior_MCMC, 
                                                                 y = y, sigma = sigma)
        
        w = calc_weight(posterior_MCMC_sample = posterior_MCMC$posterior, 
                        posterior_target_sample = log_posterior_target_sample$posterior, 
                        mu_posterior_MCMC_sample = posterior_MCMC$post[,'mu'])
        
        IS_Est_mu[i,j] = w$IS_Est_mu
        ESS_IS[i,j] = w$ESS_IS
        ESS_MCMC[i,j] = w$ESS_MCMC
        ESS[i,j] = w$ESS
      }
    }
    output = list(hyp_par = posterior_MCMC$hyp_par, mean_mu_posterior_MCMC_sample = mean(posterior_MCMC$post[,'mu']), 
                  IS_Est_mu = IS_Est_mu, ESS_IS = ESS_IS,  ESS_MCMC = ESS_MCMC, ESS = ESS)
    return(output = output)
}


######################################################################################################
# Example
L = c(length(mu_interval),length(tau_interval))
#L = c(141,33)
target_ESS = 5000


par_stan2 = c(5,7)

MCMC_stan2 = get_log_posterior_MCMC_sample(mu0 = par_stan2[1], tau0 = par_stan2[2], y = y, sigma = sigma, n_iter_stan = n_iter_stan, target_ESS = 5000)

output2 = all_target_sample_weight_given_a_MCMC_sample(low_mu0 = low_mu0, high_mu0 = high_mu0,
                                               low_tau0 = low_tau0, high_tau0 = high_tau0,
                                               L = L, y = y, sigma = sigma, posterior_MCMC = MCMC_stan2)

save(output2, file = "output2.Rdata")
load("output2.Rdata")


zlim = c(0,range(t(output2$ESS), output2$ESS_IS, t(output2$ESS_MCMC))[2])
cols <- cm.colors(3)
below.colors <- colorRampPalette(c(cols[1],cols[2]))
above.colors <- colorRampPalette(c(cols[2],cols[3]))
zz <- pretty(zlim,50)
cols=c(below.colors(sum(zz<=target_ESS)),above.colors(sum(zz>target_ESS)))

### Figure ESS
par(mar = c(5, 5, 3, 3))
filled.contour(x = mu_interval, y = tau_interval, z = t(output2$ESS), 
               zlim = zlim,
               plot.title = title(main = "ESS", xlab = expression(mu[0]), ylab = expression(tau[0]), cex.main = 2.25, cex.lab = 2), 
               levels = zz, 
               col = cols,
               plot.axes={axis(1, cex.axis = 2); axis(2, cex.axis = 2); points(par_stan2[1], par_stan2[2])},
               key.axes = axis(4, cex.axis = 1.5),
               par.settings = list(axis.text=list(fontfamily = "LM Roman 10"),
                                 par.xlab.text=list(fontfamily = "LM Roman 10"),
                                 par.ylab.text=list(fontfamily = "LM Roman 10"),
                                 par.main.text=list(fontfamily = "LM Roman 10"),                                                                             
                                 par.sub.text=list(fontfamily = "LM Roman 10")
               ))
               
### Figure ESS_MCMC
par(mar = c(5, 5, 3, 3))
filled.contour(x = mu_interval, y = tau_interval, z = t(output2$ESS_MCMC), 
               zlim = zlim,
               plot.title = title(main = "ESS_MCMC", xlab = expression(mu[0]), ylab = expression(tau[0]), cex.main = 2.25, cex.lab = 2), 
               levels = zz,
               col = cols,
               plot.axes={axis(1, cex.axis = 2); axis(2, cex.axis = 2); points(par_stan2[1], par_stan2[2])},
               key.axes = axis(4, cex.axis = 1.5),
               par.settings = list(axis.text=list(fontfamily = "LM Roman 10"),
                                   par.xlab.text=list(fontfamily = "LM Roman 10"),
                                   par.ylab.text=list(fontfamily = "LM Roman 10"),
                                   par.main.text=list(fontfamily = "LM Roman 10"),                                                                             
                                   par.sub.text=list(fontfamily = "LM Roman 10")
               ))


### Figure ESS_IS
par(mar = c(5, 5, 3, 3))
filled.contour(x = mu_interval, y = tau_interval, z = t(output2$ESS_IS),
               zlim = zlim,
               plot.title = title(main = "ESS_IS", xlab = expression(mu[0]), ylab = expression(tau[0]), cex.main = 2.25, cex.lab = 2), 
               levels = zz,
               col = cols,
               plot.axes={axis(1, cex.axis = 2); axis(2, cex.axis = 2); points(par_stan2[1], par_stan2[2])},
               key.axes = axis(4, cex.axis = 1.5),
               par.settings = list(axis.text=list(fontfamily = "LM Roman 10"),
                                   par.xlab.text=list(fontfamily = "LM Roman 10"),
                                   par.ylab.text=list(fontfamily = "LM Roman 10"),
                                   par.main.text=list(fontfamily = "LM Roman 10"),                                                                             
                                   par.sub.text=list(fontfamily = "LM Roman 10")
               ))


