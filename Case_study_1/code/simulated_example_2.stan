data {
  int<lower=0> N;                   // sample size
  // int<lower=0> N_t;
  // int<lower=0> N_c;
  // int<lower=0, upper=N> ii_t[N_t];
  // int<lower=0, upper=N> ii_c[N_c];
  vector[N] y;                      // observed outcome
  vector[N] w;                      // treatment assigned
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
   y ~ normal(alpha + tau * w, sigma_t*w + sigma_c*(1 - w));
}
generated quantities{
  
  real tau_fs;                      // finite-sample average treatment effect  
  
  vector[2] sci_tab;                // science table
  vector[2] mu;                     // location vector mu
  matrix[2, 2] Sigma;               // covariance matrix Sigma
  
  mu[1] = alpha;
  mu[2] = alpha + tau;
  
  Sigma[1, 1] = sigma_c^2;
  Sigma[1, 2] = rho*sigma_c*sigma_t;
  Sigma[2, 1] = rho*sigma_c*sigma_t;
  Sigma[2, 2] = sigma_t^2;
  
  sci_tab = multi_normal_rng(mu, Sigma);
  
  // sci_tab[ii_c, 1] = y[ii_c];
  // sci_tab[ii_t, 2] = y[ii_t];
  tau_fs = sci_tab[2] - sci_tab[1];
}

