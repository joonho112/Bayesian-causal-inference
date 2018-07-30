data {
  int<lower=0> N;                   // sample size
  real y[N];                        // observed outcome
  int<lower=0,upper=1> w[N];        // treatment assigned
  real<lower=0,upper=1> rho;        // assumed correlation between the potential outcomes
}
parameters {
  real alpha;                       // intercept
  real tau;                         // super-population average treatment effect
  real<lower=0, upper=10> sigma_c;  // residual SD for the control
  real<lower=0, upper=10> sigma_t;  // residual SD for the treated
}
model {
   // PRIORS
   alpha ~ normal(0, 5);            
   tau ~ normal(0, 5);
   sigma_c ~ cauchy(0, 5);          
   sigma_t ~ cauchy(0, 5);

   // LIKELIHOOD
   for(n in 1:N) 
       y[n] ~ normal(alpha + tau * w[n], sigma_t*w[n] + sigma_c*(1 - w[n]));
}
generated quantities{
  real tau_fs;                      // finite-sample average treatment effect  
  real y0[N];                       // potential outcome if W = 0
  real y1[N];                       // potential outcome if W = 1
  real tau_unit[N];                 // unit-level treatment effect
  for(n in 1:N){
    real mu_c = alpha;            
    real mu_t = alpha + tau;      
    if(w[n] == 1){                
      y0[n] = mu_c + (sigma_c/sigma_t)*(y[n] - mu_t); 
      y1[n] = y[n];
    }else{                        
      y0[n] = y[n];       
      y1[n] = mu_t + (sigma_t/sigma_c)*(y[n] - mu_c); 
    }
    tau_unit[n] = y1[n] - y0[n];
  }
  tau_fs = mean(tau_unit);        
}

