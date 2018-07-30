functions {  // Define quantile function (inverse CDF)
  
  real quantile(vector x, real p){
    
    int n;            // the length of vector x
    real index;       // the integer index of p
    int lo;           // the lower integer cap of index
    int hi;           // the higher integer cap of index
    real h;           // the generated weight between lo and hi
    real qs;          // the weighted average of x[lo] and x[hi]
    
    // set index, floor, and ceiling
    n = num_elements(x);           
    index = 1 + (n - 1)*p;             
    
    lo = 1; 
    while ((lo + 1) < index) 
      lo = lo + 1; 
      
    hi = lo + 1; 

    // calculate quantile
    h = index - lo; 
    qs = (1 - h)*sort_asc(x)[lo] + h*sort_asc(x)[hi];     
    
    return qs;                
  }
}
data {
  int<lower=0> N;                  // sample size
  int<lower=0> N_cov;              // number of covariates

  real y[N];                       // continuous outcome
  int<lower=0,upper=1> z[N];       // treatment assigned
  vector[N_cov] x[N];              // covariates
  vector[N_cov] xz_inter[N];       // interaction terms
  
  real<lower=0,upper=1> rho;       // assumed correlation between the potential outcomes
}

parameters {
  real alpha;                       // intercept
  row_vector[N_cov] beta;           // coefficients for x[N]
  row_vector[N_cov] beta_inter;     // coefficients for x_inter[N] 
  
  real tau;                         // treatment effect
  real<lower=0, upper=10> sigma_t;  // residual SD for the treated
  real<lower=0, upper=10> sigma_c;  // residual SD for the control
}

model {
   
   // PRIORS
   alpha ~ normal(0, 100);           // unscaled normal priors for coefficients
   beta ~ normal(0, 100);
   beta_inter ~ normal(0, 100);
   tau ~ normal(0, 100);
   sigma_t ~ inv_gamma(1, 0.01);     // inverse gamma priors for scale parameters
   sigma_c ~ inv_gamma(1, 0.01);

   // LIKELIHOOD
   for(n in 1:N)
     y[n] ~ normal(alpha + beta * x[n] + beta_inter * xz_inter[n] + tau * z[n], sigma_t*z[n] + sigma_c*(1-z[n]));
}

// FOR FINITE SAMPLE INFERENCE

generated quantities{
  
  // finite sample average treatment effect (ATE) & mean values
  real y0_bar;
  real y1_bar;
  real tau_samp;   
  
  // finite sample quantile treatment effects (QTEs)
  real tau_qte25;  
  real tau_qte50;
  real tau_qte75;

  { // to create temporary variables

    // Generate science table
    real y0[N];
    real y1[N];
    real tau_ind[N];  // individual treatment effects (ITE)

    for(n in 1:N){

      // Predicted mean 
      real mu_t = alpha + beta*x[n] + beta_inter*x[n] + tau;
      real mu_c = alpha + beta*x[n];          
      
      if(z[n] == 1){    // impute Y_mis (z = 1)
        y0[n] = normal_rng(mu_c + rho*(sigma_c/sigma_t)*(y[n] - mu_t), sigma_c*sqrt(1 - rho^2));
        y1[n] = y[n];
      }else{            // impute Y_mis (z = 0)
        y0[n] = y[n];
        y1[n] = normal_rng(mu_t + rho*(sigma_t/sigma_c)*(y[n] - mu_c), sigma_t*sqrt(1 - rho^2));
      }
      
      tau_ind[n] = y1[n] - y0[n];  // generate ITE   
        
    } // end for loop

    // store mean values and ATE
    y0_bar = mean(y0);
    y1_bar = mean(y1);
    tau_samp = mean(tau_ind);
    
    // store QTEs
    tau_qte25 = quantile(to_vector(y1), 0.25) - quantile(to_vector(y0), 0.25); 
    tau_qte50 = quantile(to_vector(y1), 0.50) - quantile(to_vector(y0), 0.50); 
    tau_qte75 = quantile(to_vector(y1), 0.75) - quantile(to_vector(y0), 0.75); 
    

  } // end temporary variables
}
