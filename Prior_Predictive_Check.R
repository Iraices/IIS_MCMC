# Prior Predictive Check

## load data
source('Load_biomanipulation_data.R')
## generate candidate priors

## typical kind of expert judgement one might have
y_int <- c(-10, 20)

# y_int1 <- c(10, 40)

# threshold 
th_1 <- 0.70

filterout_priors <- function(low_mu0 = -100,high_mu0 = 100,
                      low_tau0 = 5, high_tau0 = 50, y_int = y_int, length.out = 100,
                      th_1 = 0.70){
val <- expand.grid(mu0  = seq(from = low_mu0, to = high_mu0, length.out = length.out),#by = 0.25),
                   tau0 = seq(from = low_tau0, to = high_tau0, length.out = length.out))# by = 0.25))
## sample data for each choice of prior
val$p <- unlist(lapply(1:nrow(val),function(i){
  n_iter <- 10000 
 
  tau_mu <- runif(n_iter, min = 1, max = val$tau0[i])
  mu <- rnorm(n_iter, mean = val$mu0[i], sd = tau_mu)
  k <- runif(n_iter, min = 1, max = 5)
  tau_study <- tau_mu * k
  y_gen <- rnorm(n_iter, mean = mu, sd = sqrt(tau_study^2))
  
  mean(y_gen > y_int[1] & y_gen < y_int[2])
}))
## select sets of priors
#th_1 = 0.75
val$dichotomic <- val$p > th_1
## Select set of priors to use mu0 and tau0
val_dichotomic_TRUE <- filter(val, dichotomic == TRUE)
return(list(mu0_set_prior = range(val_dichotomic_TRUE$mu0),
       tau0_set_prior = range(val_dichotomic_TRUE$tau0),
       val=val))
}

out <- filterout_priors(low_mu0 = -30,high_mu0 = 30,
                        low_tau0 = 5,high_tau0 = 50, 
                        y_int = y_int,length.out = 100, th_1 = 0.70)
levelplot(p ~ mu0*tau0, data = out$val, xlab = 'mu0', ylab = 'tau0')
levelplot(dichotomic ~ mu0*tau0, data = out$val, xlab = 'mu0', ylab = 'tau0')
out$mu0_set_prior
out$tau0_set_prior


out <- filterout_priors(low_mu0 = -6,
                        high_mu0 = 16,
                        low_tau0 = 5,
                        high_tau0 = 10, 
                        y_int = y_int,
                        length.out = 200,th_1 = 0.70)


levelplot(p ~ mu0*tau0, data = out$val, xlab = 'mu0', ylab = 'tau0',
          par.settings=list(axis.text=list(fontfamily="LM Roman 10"),
                            par.xlab.text=list(fontfamily="LM Roman 10"),
                            par.ylab.text=list(fontfamily="LM Roman 10"),
                            par.main.text=list(fontfamily="LM Roman 10"),                                                                             
                            par.sub.text=list(fontfamily="LM Roman 10")
          ))   


plot_1 <- levelplot(dichotomic ~ mu0*tau0, data = out$val, cut = 2, breaks = 0.70, at = c(0,0.70,1), 
                   scales = list(tck = c(1,0), x = list(cex = 2), y = list(cex = 2)),
                   xlab = list(label = expression(mu[0]), cex = 2) , ylab = list(label = expression(tau[0]), cex = 2),
                   colorkey = list(labels = list(cex = 1.75)),
                   par.settings = list(axis.text=list(fontfamily = "LM Roman 10"),
                                     par.xlab.text=list(fontfamily = "LM Roman 10"),
                                     par.ylab.text=list(fontfamily = "LM Roman 10"),
                                     par.main.text=list(fontfamily = "LM Roman 10"),                                                                             
                                     par.sub.text=list(fontfamily = "LM Roman 10")
                   ))
                   
out$mu0_set_prior
out$tau0_set_prior


## range for sets of priors to use

mu0_set_prior <- round(out$mu0_set_prior)
tau0_set_prior <- round(out$tau0_set_prior)
save(mu0_set_prior,tau0_set_prior,y_int,
     file='setofpriorstouse.Rdata')

########################################################################################################################

## plot the distribution for these priors compared to the data
val <- out$val
i = which(val$p > th_1)[1]
val[i,]

n_iter <- 10000

## Figure Prior Predictive Check 
par(mar = c(5, 5, 3, 3), family = "LM Roman 10")
plot(x = c(-200,200), y = c(0,1), type='n', ylab = 'cdf', xlab = 'Response variable', cex.lab = 2, cex.axis = 2)

 
 red_val <- val[val$p > th_1,]
 
 for(i in sample.int(nrow(red_val),50)){
  tau_mu <- runif(n_iter, min = 1, max = red_val$tau0[i])
  mu <- rnorm(n_iter, mean = red_val$mu0[i], sd = tau_mu)
  k <- runif(n_iter, min = 1, max = 5)
  tau_study <- tau_mu * k
  y_gen <- rnorm(n_iter, mean = mu, sd = sqrt(tau_study^2))
  
  lines(sort(y_gen),(1:length(y_gen))/length(y_gen),col='red')
 }
 
segments(x0 = y_int[1], x1 = y_int[2], y0 = 0, col='black',lwd=2) 
segments(x0 = y_int[1],  y0 = 0, y1 = 1, col='black',lwd=2) 
segments(x0 = y_int[2],  y0 = 0, y1 = 1, col='black',lwd=2) 
segments(x0 = y_int[1], x1 = y_int[2], y0 = 1, col='black',lwd=2) 


#################################################################

border_values <- data.frame(x=out$val$mu0[abs(out$val$p-0.70)<=0.01],y=out$val$tau0[abs(out$val$p-0.70)<=0.01])

border_values_edge_df <- border_values %>% 
  mutate(bin_1=cut(x, 100, labels = FALSE)) %>% 
  group_by(bin_1) %>% 
  mutate(q=rank(y)/n()) %>% 
  filter(q>0.98) %>% 
  ungroup() 

fit_model <- lm(y ~ poly(x, 2, raw = TRUE), data = border_values_edge_df)

data_order <- data.frame(x = border_values_edge_df$x, y= predict(fit_model, border_values_edge_df))
data_order <- data_order[order(data_order$x),]

plot_2 <- xyplot(predict(fit_model, data_order) ~ data_order$x, col = 'black', type = 'l', xlim = c(-6,16), ylim = c(4,10), 
                xlab = expression(mu[0]), ylab = expression(tau[0]))

## Figure levelplot
plot_1 + as.layer(plot_2)


########################################################################################

plot(x = border_values_edge_df$x, y = predict(fit_model, border_values_edge_df), col = 'black',  xlim = c(-6,16), ylim = c(5,10), xlab = 'mu0', ylab = 'tau0')
par(new = TRUE)
plot(border_values_edge_df$x, border_values_edge_df$y, xlim = c(-6,16), ylim = c(5,10), xlab = 'mu0', ylab = 'tau0')


###################################

change_output <- tidy(fit_model)
 
 region <- function(x, coef){
   l <- length(coef)
   grade <- seq(from = 0,to = l, by = 1)
   func <- 0
   for (i in 1:l){
     func <- func + coef[i] * x^(grade[i]) 
   }
   return(func)
 }
 
 reg <- function(x){region(x, coef =  change_output$estimate)}
 
 
 