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
  
  real y[N];                       // continuous outcome
  int<lower=0,upper=1> z[N];       // treatment assigned
  
  real<lower=0,upper=1> rho;       // assumed correlation between the potential outcomes
}

parameters {
  real alpha;                       // intercept
  real tau;                         // treatment effect
  real<lower=0, upper=10> sigma;    // residual SD for the treated and the control
}

model {
   
   // PRIORS
   alpha ~ normal(0, 100);           // unscaled normal priors for coefficients
   tau ~ normal(0, 100);
   sigma ~ inv_gamma(1, 0.01);       // inverse gamma for scale parameter

   // LIKELIHOOD
   for(n in 1:N)
     y[n] ~ normal(alpha + tau * z[n], sigma);
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

      real mu_t = alpha + tau;      // mean for the treated
      real mu_c = alpha;            // mean for the control
      
      if(z[n] == 1){    // impute Y_mis (z = 1)
        y0[n] = mu_c + rho*(y[n] - mu_t);
        y1[n] = y[n];
      }else{            // impute Y_mis (z = 0)
        y0[n] = y[n];
        y1[n] = mu_t + rho*(y[n] - mu_c);
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

