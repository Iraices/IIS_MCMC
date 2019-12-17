/// Stan code 
data {
  int<lower=2> N;          // Number of studies
  real y[N];               // Observed intervention effect 
  real<lower=0> sigma[N];  // Within-study variance (SD)
  real mu0;		   // hyperparameter for the prior on mu
  real<lower=0>	tau0;	   // hyperparameter for the prior on tau 
 }

transformed data {
  real<lower=0> sigma2[N];
  for (i in 1:N){
    sigma2[i] = sigma[i] * sigma[i];
  }
}

parameters {                                
  real mu;                   // Summary intervention effect   
  real<lower=0> tau_mu;      // SD of summary intervention effect    
  real log_k;                // Log-proportional constant between SD of Summary intervention effect  and between-study variation 
}

transformed parameters {
  real<lower=0> tau_overall[N];  // Sum of Within-study variance (SD) and between-study variance (SD)
  real<lower=0> tau_study;       // Between-study variance (SD)
  real<lower=0> k;		 // Log-proportional constant between SD of Summary intervention effect and between-study variation 
  k = exp(log_k); 
  tau_study = tau_mu * k;
  for (i in 1:N){
    tau_overall[i] = sqrt(sigma2[i] + tau_study * tau_study);
   }
}

model {
 tau_mu ~ normal(0,tau0); 
 mu ~ normal(mu0, tau_mu);
 log_k ~ normal(0,1.5);
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
 log_posterior = log_lik + normal_lpdf(mu| mu0, tau_mu) + normal_lpdf(tau_mu| 0, tau0) + normal_lpdf(log_k | 0, 1.5);
}
