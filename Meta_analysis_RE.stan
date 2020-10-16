/// Stan code 
data {
  int<lower=2> N;          // Number of studies
  real y[N];               // Observed intervention effect 
  real<lower=0> sigma[N];  // Within-study variance (SD)
  real mu0;		   // hyperparameter for the prior on mu
  real<lower=2>	tau0;	   // hyperparameter for the prior on tau 
 }

transformed data {
  real<lower=0> sigma2[N];
  for (i in 1:N){
    sigma2[i] = sigma[i] * sigma[i];
  }
}

parameters {                                
  real mu;                   // Summary intervention effect   
  real<lower = 1, upper = tau0> tau_mu;      // SD of summary intervention effect    
  real<lower = 1, upper = 5> k;           // proportional constant between SD of Summary intervention effect  and between-study variation 
}

transformed parameters {
  real<lower=0> tau_overall[N];  // Sum of Within-study variance (SD) and between-study variance (SD)
  real<lower=1> tau_study;       // Between-study variance (SD) 
  tau_study = tau_mu * k;
  for (i in 1:N){
    tau_overall[i] = sqrt(sigma2[i] + tau_study * tau_study);
   }
}

model {
 tau_mu ~ uniform(1, tau0); 
 mu ~ normal(mu0, tau_mu);
 k ~ uniform(1,5);
  for (i in 1:N){
   	y[i] ~ normal(mu, tau_overall[i]);
  }
}

generated quantities {
 real log_lik; 		 // log_likelihood
 real log_posterior;
 vector[N] log_lik_data; 
  for (i in 1:N){
    log_lik_data[i] = normal_lpdf(y[i]|mu, tau_overall[i]);
  }
 log_lik = sum(log_lik_data);
 log_posterior = log_lik + normal_lpdf(mu| mu0, tau_mu) + uniform_lpdf(tau_mu| 1, tau0) + uniform_lpdf(k| 1, 5);
}

