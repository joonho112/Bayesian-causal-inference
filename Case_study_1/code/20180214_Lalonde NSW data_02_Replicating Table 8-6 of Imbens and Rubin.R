
###'######################################################################
###'
###' Bayesian Model-Based Inference for Completely Randomized Experiments
###' 
###' 
###' Replicating Table 8.6 of Imbens & Rubin's Textbook 
###' 
###' 
###' 20180228 JoonHo Lee
###'
###'
###'


###'######################################################################
###'
###' Basic settings
###'
###'


### Remove previous workspace

rm(list=ls())


### Set working directory 

setwd("~/Stan/20180204_Bayesian Model-Based Inference for RCT")


### Call libraries

library(readxl)
library(foreign)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

library(rstan)
library(shinystan)
library(Matching)


# Set Stan options
options(mc.cores = parallel::detectCores())




###'######################################################################
###'
###' Import dataset & summary statistics
###'
###'

### Use example dataset form Matching package
data(lalonde, package = "Matching")


### Summary statistics
sumstats <- lalonde %>%
  
  # Find the mean, std. dev., min, and max for each variable 
  summarise_all(funs(mean, sd, min, max)) %>%
  
  # Move summary stats to columns
  gather(key, value, everything()) %>% 
  separate(key, into = c("variable", "stat"), sep = "_") %>%
  spread(stat, value) %>%
  dplyr::select(variable, mean, sd, min, max) %>%
  
  # Round all numeric variables to two decimal point
  mutate_each(funs(round(., 2)), -variable)

knitr::kable(sumstats, caption = "Summary statistics")




# ###'######################################################################
# ###'
# ###' Setting up dataset: #Option 1. Standardizing to improve MCMC sampling
# ###'
# ###'
# 
# ### Standardize age covariate
# lalonde$age_std <- (lalonde$age - mean(lalonde$age))/sd(lalonde$age)
# 
# ### Log-transform all earnings variables
# lalonde$log_re74 <- log(lalonde$re74 + 1)
# lalonde$log_re75 <- log(lalonde$re75 + 1)
# lalonde$log_re78 <- log(lalonde$re78 + 1)
# 
# ### Set outcome variable & standardize to improve MCMC sampling
# y <- lalonde$log_re78
# y <- (y - mean(y))/sd(y)
# 
# ### Set treatment variable
# z <- lalonde$treat
# 
# ### Collect coviriates x
# x <- as.matrix(lalonde[, c("age_std", "nodegr", "black", "hisp", "married", 
#                            "log_re74", "u74", "log_re75", "u75")])





###'######################################################################
###'
###' Setting up dataset: #Option 2. Earnings in thousands
###' 
###' To replicate Table 8.6 in Imbens & Rubin's Chapter 8
###'
###'

### Change the scales of all earnings variables in thousands
lalonde$re74 <- lalonde$re74/1000
lalonde$re75 <- lalonde$re75/1000
lalonde$re78 <- lalonde$re78/1000

### Set outcome variable & standardize to improve MCMC sampling
y <- lalonde$re78

### Set treatment variable
z <- lalonde$treat

### Collect coviriates x
x <- as.matrix(lalonde[, c("age", "educ", "married", "nodegr", "black", 
                           "re74", "u74", "re75", "u75")])


###'######################################################################
###'
###' Neyman's approach
###' 
###' : Simple difference in means
###' 1.794 (0.632)
###'
###'

summary(lm(y ~ z))




###'######################################################################
###'
###' Model 1. General Form without Covariates_Perfect correlation_Common variance
###' 
###' (1) Mean Covariate Dependent: No
###' 
###' (2) Variance Treatment Specific: No
###' 
###' (3) Potential Outcome Independent: No
###' 
###' (4) Two-Part Model: No     
###'
###'


### Set correlation between potential outcomes

rho_value <- 1.0



### Assign data list for stan

stan_data <- list(  y = y,
                    z = z,
                    x = x,
                    N = nrow(lalonde),
                    N_cov = ncol(x), 
                    rho = rho_value )


### Fit stan model

stanfit_mod1 <- stan(file = "code/Model_05_General Form without Covariates_Perfect correlation_Common residual variance.stan",
                     data = stan_data,
                     iter = 1000, chains = 4)

print(stanfit_mod1)





###'######################################################################
###'
###' Model 2. General Form without Covariates
###' 
###' (1) Mean Covariate Dependent: No  => Without covariates
###' 
###' (2) Variance Treatment Specific: Yes  => heterogeneous variances 
###' 
###' (3) Potential Outcome Independent: Yes  => rho = 0.0 
###' 
###' (4) Two-Part Model: No     
###'
###'


### Set correlation between potential outcomes

rho_value <- 0.0



### Assign data list for stan

stan_data <- list(  y = y,
                    z = z,
                    x = x,
                    N = nrow(lalonde),
                    N_cov = ncol(x), 
                    rho = rho_value )


### Fit stan model

stanfit_mod2 <- stan(file = "code/Model_01_General Form without Covariates.stan",
                     data = stan_data,
                     iter = 1000, chains = 4)

print(stanfit_mod2)






###'######################################################################
###'
###' Model 3. General Form with Covariates
###' 
###' (1) Mean Covariate Dependent: Yes  => With covariates
###' 
###' (2) Variance Treatment Specific: Yes  => heterogeneous variances 
###' 
###' (3) Potential Outcome Independent: Yes  => rho = 0.0 
###' 
###' (4) Two-Part Model: No     
###'
###'


### Set correlation between potential outcomes

rho_value <- 0.0



### Add mean-centered covariates and their interaction terms 

x_mean <- apply(x, 2, mean, na.rm = TRUE)

x_mean_mat <- matrix(rep(x_mean, each = nrow(x)), nrow = nrow(x))

x_c_mat <- x - x_mean_mat  # mean-centering

xz_inter <- x_c_mat*z         # interaction terms (elementwise multiplication)

colnames(xz_inter) <-  paste0(colnames(x), "_z")



### Assign data list for stan

stan_data <- list(  y = y,
                    z = z,
                    x = x_c_mat,
                    xz_inter = xz_inter,  
                    N = nrow(lalonde),
                    N_cov = ncol(x_c_mat), 
                    rho = rho_value )


### Fit stan model

stanfit_mod3 <- stan(file = "code/Model_02_General Form with Covariates.stan",
                     data = stan_data,
                     iter = 1000, chains = 4)

print(stanfit_mod3)






###'######################################################################
###'
###' Model 4. Two-Part Model with Covariates
###' 
###' (1) Mean Covariate Dependent: Yes  => With covariates
###' 
###' (2) Variance Treatment Specific: Yes  => heterogeneous variances 
###' 
###' (3) Potential Outcome Independent: Yes  => rho = 0.0 
###' 
###' (4) Two-Part Model: Yes => Two-Part model     
###'
###'


### Set correlation between potential outcomes

rho_value <- 0.0



### Add mean-centered covariates and their interaction terms 

x_mean <- apply(x, 2, mean, na.rm = TRUE)

x_mean_mat <- matrix(rep(x_mean, each = nrow(x)), nrow = nrow(x))

x_c_mat <- x - x_mean_mat  # mean-centering

xz_inter <- x_c_mat*z         # interaction terms (elementwise multiplication)

colnames(xz_inter) <-  paste0(colnames(x), "_z")



### Assign data list for stan

stan_data <- list(  y = y,
                    z = z,
                    x = x_c_mat,
                    xz_inter = xz_inter,  
                    y_pos = as.numeric(lalonde$re78 > 0), ## for two-part model
                    N = nrow(lalonde),
                    N_cov = ncol(x_c_mat), 
                    rho = rho_value )


### Fit stan model

stanfit_mod4 <- stan(file = "code/Model_06_General Form with Covariates_Two-Part Model.stan",
                     data = stan_data,
                     iter = 1000, chains = 4)

print(stanfit_mod4)




###'######################################################################
###'
###' Save estimation results
###'
###'

stanfit_collect <- list(stanfit_mod1, 
                        stanfit_mod2, 
                        stanfit_mod3, 
                        stanfit_mod4)

save(stanfit_collect, file = "data/stanfit_collect_Table8-6.RData")




###'######################################################################
###'
###' Generate Table 8.6
###'
###'

### Load the saved RData object

load(file = "data/stanfit_collect_Table8-6.RData")



### Generate empty matrix to collect estimates

tau_est_collect <- matrix(ncol = 4, nrow = 0)



### Collect tau estimates

for (i in 1:4){
  
  tau_est <- summary(stanfit_collect[[i]], 
                     c("tau_samp", "tau_qte25", "tau_qte50", "tau_qte75"))$summary
  
  
  tau_mean <- sprintf("%.2f", round(as.vector(tau_est[, 1]), 2)) 
  tau_sd <- paste0("(", sprintf("%.2f", round(as.vector(tau_est[, 3]), 2)), ")")
  
  
  est <- cbind.data.frame(tau_mean, tau_sd) %>% 
    unite(est, tau_mean, tau_sd, sep = " ")
  
  tau_est_collect <- rbind(tau_est_collect, est$est)

}


### Generate Table 8.6

mod_info_vec <- c(FALSE, FALSE, FALSE, FALSE, 
                  FALSE, TRUE, TRUE, FALSE, 
                  TRUE, TRUE, TRUE, FALSE, 
                  TRUE, TRUE, TRUE, TRUE)

mod_info_mat <- matrix(mod_info_vec, ncol = 4, nrow = 4, byrow = TRUE)

table8_6 <- cbind.data.frame(mod_info_mat, tau_est_collect)

names(table8_6) <- c("Covariates", "Specific_Variance", "Independent_PO", "Two_Part", 
                     "ATE_Sample", "QTE25%", "QTE50%", "QTE75%")


### Save the resulting table

write.csv(table8_6, file = "table/table8.6 replication.csv")






