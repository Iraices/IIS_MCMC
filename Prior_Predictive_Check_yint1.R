# Prior Predictive Check

## load data
source('Load_biomanipulation_data.R')
## generate candidate priors

## typical kind of expert judgement one might have
y_int = c(10, 40)

# threshold 
th_1 = 0.75

filterout_priors <- function(low_mu0 = -100,high_mu0 = 100,
                      low_tau0 = 1,high_tau0 = 20, y_int = y_int, length.out = 10,
                      th_1 = 0.75){
val = expand.grid(mu0  = seq(from = low_mu0, to = high_mu0, length.out = length.out),#by = 0.25),
                   tau0 = seq(from = low_tau0, to = high_tau0,length.out = length.out))# by = 0.25))
## sample data for each choice of prior
val$p = unlist(lapply(1:nrow(val),function(i){
  n_iter = 10000 
 
  tau_mu = abs(rnorm(n_iter, mean = 0, sd = val$tau0[i]))
  mu = rnorm(n_iter, mean = val$mu0[i], sd = tau_mu)
  k = exp(rnorm(n_iter, mean = 0, sd = 1.5))
  tau_study = tau_mu * k
  y_gen = rnorm(n_iter, mean = mu, sd = sqrt(tau_study^2))
  
  mean(y_gen > y_int[1] & y_gen < y_int[2])
}))
## select sets of priors
val$dichotomic <- val$p > th_1
## Select set of priors to use mu0 and tau0
val_dichotomic_TRUE = filter(val, dichotomic == TRUE)
return(list(mu0_set_prior = range(val_dichotomic_TRUE$mu0),
       tau0_set_prior = range(val_dichotomic_TRUE$tau0),
       val=val))
}

out <- filterout_priors(low_mu0 = 0,high_mu0 = 50,
                        low_tau0 = 1,high_tau0 = 9, 
                        y_int = y_int,length.out = 100, th_1=0.75)
levelplot(p ~ mu0*tau0, data = out$val, xlab = 'mu0', ylab = 'tau0')
levelplot(dichotomic ~ mu0*tau0, data = out$val, xlab = 'mu0', ylab = 'tau0')
out$mu0_set_prior
out$tau0_set_prior


out <- filterout_priors(low_mu0 = 10,high_mu0 = 40,
                        low_tau0 = 1,high_tau0 = 9, 
                        y_int = y_int,length.out = 100, th_1=0.75)
levelplot(p ~ mu0*tau0, data = out$val, xlab = 'mu0', ylab = 'tau0')
levelplot(dichotomic ~ mu0*tau0, data = out$val, xlab = 'mu0', ylab = 'tau0')
out$mu0_set_prior
out$tau0_set_prior

plot_1 = levelplot(dichotomic ~ mu0*tau0, data = out$val, cut = 2, breaks = 0.75, at = c(0,0.75,1), xlab = 'mu0', ylab = 'tau0')

## range for sets of priors to use

mu0_set_prior = round(out$mu0_set_prior)
tau0_set_prior = round(out$tau0_set_prior)
save(mu0_set_prior,tau0_set_prior,y_int,
     file='setofpriorstouse1.Rdata')

########################################################################################################################

## plot the distribution for these priors compared to the data
val <- out$val
i = which(val$p > th_1)[1]
val[i,]

n_iter = 10000

png('CDF111.png')

 plot(c(-200,200),c(0,1),type='n', ylab = 'cdf', xlab = 'Response variable')
 
 red_val <- val[val$p > th_1,]
 
 for(i in sample.int(nrow(red_val),50)){
  tau_mu = abs(rnorm(n_iter, mean = 0, sd = red_val$tau0[i]))
  mu = rnorm(n_iter, mean = red_val$mu0[i], sd = tau_mu)
  k = exp(rnorm(n_iter, mean = 0, sd = 1.5))
  tau_study = tau_mu * k
  y_gen = rnorm(n_iter, mean = mu, sd = sqrt(tau_study^2))
  
  lines(sort(y_gen),(1:length(y_gen))/length(y_gen),col='red')
 }
 legend('topleft',c('elicited range','generated data'),col = c('black','red'),lty=c(1,1))
 
segments(x0 = y_int[1], x1 = y_int[2], y0 = 0, col='black',lwd=2) 
segments(x0 = y_int[1],  y0 = 0, y1 = 1, col='black',lwd=2) 
segments(x0 = y_int[2],  y0 = 0, y1 = 1, col='black',lwd=2) 
segments(x0 = y_int[1], x1 = y_int[2], y0 = 1, col='black',lwd=2) 

dev.off()


#################################################################

border_values = data.frame(x=out$val$mu0[abs(out$val$p-0.75)<=0.01],y=out$val$tau0[abs(out$val$p-0.75)<=0.01])

border_values_edge_df <- border_values %>% 
  mutate(bin_1=cut(x, 100, labels = FALSE)) %>% 
  group_by(bin_1) %>% 
  mutate(q=rank(y)/n()) %>% 
  filter(q>0.98) %>% 
  ungroup() 

fit_model <- lm(y ~ poly(x, 2, raw = TRUE), data = border_values_edge_df)

change_output = tidy(fit_model)
 
 region = function(x, coef){
   l = length(coef)
   grade = seq(from = 0,to = l, by = 1)
   func = 0
   for (i in 1:l){
     func = func + coef[i] * x^(grade[i]) 
   }
   return(func)
 }
 
 reg = function(x){region(x, coef =  change_output$estimate)}
 
########################################################################
load('setofpriorstouse1.Rdata')

 x = seq(range(mu0_set_prior)[1] - 1, range(mu0_set_prior)[2] + 1, 0.1)
 plot_2 = xyplot(reg(x) ~ x, col = 'black', type = 'l', xlim = c(-10,20), ylim = c(1,9), xlab = 'mu0', ylab = 'tau0')
 
 png('levelplot_yint1.png')
 plot_1 + as.layer(plot_2)
 dev.off()


 ########################################################################
  
 low_mu0 = mu0_set_prior[1]
 high_mu0 = mu0_set_prior[2]
 low_log_tau0 = log(tau0_set_prior[1]) 
 high_log_tau0 = log(tau0_set_prior[2])
 
 n_iter_stan = 20000
 mu_interval  = seq(from = low_mu0, to = high_mu0, by = 0.2)
 tau_interval = seq(from = exp(low_log_tau0), to = exp(high_log_tau0), by = 0.25)
 grid = expand.grid(mu_interval = mu_interval, tau_interval = tau_interval)
 
 ########################################################################################
 
 initial_time_yint1 = proc.time()
  
 output_orig_yint1 = estimate_lower_bound (mu0 = 25, log_tau0 = log(5), 
                                     lower_bound_mu0 = low_mu0, upper_bound_mu0 = high_mu0, 
                                     lower_bound_tau0 = low_log_tau0, upper_bound_tau0 = high_log_tau0, 
                                     y = data$y, sigma = data$sigma, y_int = y_int,
                                     n_iter_stan = n_iter_stan, region = reg, target_ESS = 5000)
 
 final_time_yint1 = proc.time() -  initial_time_yint1
 
 
 