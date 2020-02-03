data {
  int<lower=1> N;
  
  int<lower=0,upper=1> z[N];   // treatment assigned
  int<lower=0,upper=1> w[N];   // treatment received  
  int<lower=0,upper=1> y[N];   // outcomes  
}

parameters {

  // PRINCIPAL STRATUM OUTCOME MEANS
  real<lower=0,upper=1> eta_c0;
  real<lower=0,upper=1> eta_c1;

  real<lower=0,upper=1> eta_nt0;
  real<lower=0,upper=1> eta_nt1;

  
  // OVERALL PROBABILITY OF BEING A COMPLIER
  real<lower=0,upper=1> pi;

} 


model {
  
  // PRIORS FOR OUTCOME (from Imbens & Rubin (1997))
  eta_c0 ~ beta(2, 2);  
  eta_c1 ~ beta(2, 2);  

  eta_nt0 ~ beta(2, 2);  
  eta_nt1 ~ beta(2, 2);  


  // PRIORS FOR COMPLIER PROBABILITY
  // implicit prior: pi ~ Unif(0,1)


  
  // MODELS FOR OUTCOME
  for(n in 1:N){
    
    // Never Takers (treat)
    if(z[n] == 1 && w[n] == 0){
      target += log(1 - pi) + bernoulli_lpmf(y[n] | eta_nt1);
    }
    
    
    // Complier (control) or Never Taker (control)
    else if(z[n] == 0 && w[n] == 0){
      target +=  log_sum_exp(
        log(1 - pi) + bernoulli_lpmf(y[n] | eta_nt0),  // Never taker
        log(pi) + bernoulli_lpmf(y[n] | eta_c0));      // Complier (control)
    }
    
    
    // Complier (treated)
    else if(z[n] == 1 && w[n] == 1){
      target += log(pi) + bernoulli_lpmf(y[n] | eta_c1) ;  // Complier (treat)
    }
    
  }
}