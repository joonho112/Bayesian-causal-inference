functions {  // Define quantile function (inverse CDF)
  
  real quantile(vector x, real p){
    
    int n;            // the length of vector x
    real index;       // the integer index of p
    int lo;           // the lower integer cap of index
    int hi;           // the higher integer cap of index
    real h;           // the generated weight between lo and hi
    real qs;          // the weighted average of x[lo] and x[hi]
    
    // set index
    n = num_elements(x);           
    index = 1 + (n - 1)*p;             
    
    // use while loop to cast from "real" to "int""
    lo = 1; 
    while ((lo + 1) < floor(index)) 
      lo = lo + 1; 
      
    hi = 1; 
    while ((hi + 1) < ceil(index)) 
      hi = hi + 1; 

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
  real<lower=0, upper=10> sigma_t;  // residual SD for the treated
  real<lower=0, upper=10> sigma_c;  // residual SD for the control
}

model {
   
   // PRIORS
   alpha ~ normal(0, 5);           // unscaled normal priors for coefficients
   tau ~ normal(0, 5);
   sigma_t ~ cauchy(0, 5);         // half-cauchy priors for scale parameters
   sigma_c ~ cauchy(0, 5);

   // LIKELIHOOD
   for(n in 1:N)
     y[n] ~ normal(alpha + tau * z[n], sigma_t*z[n] + sigma_c*(1-z[n]));
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
        y0[n] = mu_c + rho*(sigma_c/sigma_t)*(y[n] - mu_t);
        y1[n] = y[n];
      }else{            // impute Y_mis (z = 0)
        y0[n] = y[n];
        y1[n] = mu_t + rho*(sigma_t/sigma_c)*(y[n] - mu_c);
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

