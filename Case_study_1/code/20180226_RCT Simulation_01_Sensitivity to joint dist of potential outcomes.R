
###'######################################################################
###'
###' Bayesian Model-Based Inference for Completely Randomized Experiments
###' 
###' Simulating RCT data - With covariates
###' 
###' 
###' Sensitivity analysis for 
###' the assumption of joint distribution of potential outcomes
###' 
###' 20180214 JoonHo Lee
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
###' Preparing for loop
###'  
###' Loop over rho-s
###'
###'

### Define looping variables

rho_vec <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)



### Define an array for saving stanfits

stanfit_collect <- array(list(), dim = c(length(rho_vec)))





###'######################################################################
###'
###' Simulate RCT Data With Covariates
###' 
###' Loop over 6 rho-s
###'
###'


for (i in seq(length(rho_vec))){  # (1) Loop over 6 rho-s
  
  
  ###'#######################
  ###' Assign loop variable
  ###'#######################
  
  rho_value <- rho_vec[i]
  
  
  ###'#######################
  ###' Simulate RCT Data
  ###' Assumming rho = 0.0
  ###'#######################
  
  N <- 500          # number of observations
  
  betay <- c(0.5, 1.0, -1.5)     # coefficients of X in the outcome model
  
  tau <- 0.25     # treatment effect
  
  set.seed(725)
  
  alpha <- 1.0                   # intercept
  X <- matrix(rnorm(3 * N), N)   # covariates
  
  ps <- 0.5                      # propensity score (constant in RCT)
  
  Z <- rbinom(N, 1, ps)          # treatment variable
  
  epsilon_t <- rnorm(N, 0.0, 1.0)  # error term for the treated
  epsilon_c <- rnorm(N, 0.0, 1.0)  # error term for the control
  
  Y0 <- alpha + X %*% betay + 0*tau + epsilon_c       # potential outcome (Z=0)
  Y1 <- alpha + X %*% betay + 1*tau + epsilon_t       # potential outcome (Z=1)
  Y <- Y0 * (1 - Z) + Y1 * Z                          # realization of potential outcome
  
  
  
  ###'#######################
  ###'Assign data list 
  ###'#######################
  
  stan_data <- list(y = as.vector(Y),
                    z = Z,
                    x = X,
                    N = N,                         
                    N_cov = ncol(X), 
                    rho = rho_value)
  
  
  
  ###'#######################
  ###' Fit the Stan model
  ###'#######################
  
  if (rho_value == 1.0){  # perfect correlation (rho == 1.0)
    
    stan_fit <- stan(file = "code/Model_04_General Form with Covariates_Perfect correlation.stan",
                     data = stan_data,
                     iter = 1000, chains = 4)
    
  }else{  # non-perfect correlation (rho != 1.0)
    
    stan_fit <- stan(file = "code/Model_02_General Form with Covariates.stan",
                     data = stan_data,
                     iter = 1000, chains = 4)
  }
  
  ### Store stanfit in the array 
  stanfit_collect[[i]] <- stan_fit

}    # End of loop over 6 rho-s






###'######################################################################
###'
###' Save the resulting array
###'
###'

save(stanfit_collect, file = "data/stanfit_collect_simulation.RData")




###'######################################################################
###'
###' Generate table
###' 
###' (1) Long data for plotting and tabulating
###'
###'


### Load stanfit collect

load(file = "data/stanfit_collect_simulation.RData")



###' Generate empty data frame to collect summary statistics 
###' for posterior distribution of tau

tau_summary_collect <- data.frame(matrix(ncol = 13, nrow = 0))

var_names <- c("rho", "tau", 
               "mean", "se_mean", "sd", 
               "pct2.5", "pct25", "pct50", "pct75", "pct97.5", 
               "n_eff", "Rhat")

names(tau_summary_collect) <- var_names



### Start for loop to collect summary stats

for (i in seq(length(rho_vec))){  # (1) Loop over 6 rho-s
    
  
    ### Extract stanfit object
    
    stanfit <- stanfit_collect[[i]]
    
    
    ###' Assign summary statistics for the posterior distributions of
    ###' tau and tau_samp
    
    tau_summary <- summary(stanfit, c("tau", "tau_samp", 
                                      "tau_qte25", "tau_qte50", "tau_qte75"))$summary
    
    
    ### Add labels
    
    rho_lab <- c(paste0("rho = ", sprintf("%.1f", rho_vec[i])))
    
    add_label <- matrix(rep(cbind(rho_lab), 5), nrow = 5, byrow = TRUE)
    
    tau_df <- cbind.data.frame(add_label, 
                               rownames(tau_summary), 
                               tau_summary)
    
    
    ### Assign row names and column names
    rownames(tau_df) <- c()
    
    colnames(tau_df) <- var_names
    
    
    
    ### Append to the tau_summary_collect
    
    tau_summary_collect <- rbind.data.frame(tau_summary_collect, tau_df)
}


### Check classes

sapply(tau_summary_collect, class)


### Reorder factors for tau

tau_summary_collect$tau <- factor(tau_summary_collect$tau, 
                                  levels = c("tau_samp", "tau", 
                                             "tau_qte25", "tau_qte50", "tau_qte75"), 
                                  labels = c("ATE Sample", 
                                             "ATE Population", 
                                             "QTE 25%", 
                                             "QTE 50%", 
                                             "QTE 75%"))




###'######################################################################
###'
###' Plot #1. Only ATEs: finite sample ATE vs. Super-population ATE 
###'
###'

### Prepare dataset for plotting

df_ATE <- tau_summary_collect %>%
  filter(tau %in% c("ATE Sample", "ATE Population"))



### Plot!  

p <- ggplot(data = df_ATE, aes(x = tau, group = tau, col = tau)) +
  
  # Point
  geom_point(aes(y = mean), size = 3) + 
  
  # Error bar
  geom_errorbar(aes(ymin = pct2.5, ymax = pct97.5), width = 0.2) + 
  
  
  # Value label (posterior mean & sd)
  geom_text(aes(y = mean, label = sprintf("%.3f", round(mean, 3)), 
                hjust = 0.5, vjust = -1.0)) + 
  geom_text(aes(y = mean, label = paste0("(", sprintf("%.3f", round(sd, 3)), ")"), 
                hjust = -0.5, vjust = -1.0)) + 
  
  # Faceting
  facet_grid(rho ~ .) + 
  
  # Flip coordinates
  coord_flip() +
  
  
  # Color setting
  scale_colour_manual(labels = c("Finite-sample ATE", "Super-population ATE"), 
                      values = c("firebrick1", "dodgerblue1")) + 
  
  
  # Theme settings
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank()) + 
  
  # Labels
  labs(title = "Estimated ATEs",  
       subtitle = "Posterior Means, Standard Deviations (in parentheses), and 95% Credible Intervals",  
       y = "",
       x = "Estimated Average Treatment Effect")


### Save the plot 

setwd("~/Stan/20180204_Bayesian Model-Based Inference for RCT")

ggsave("figure/Simulation_ATE_Sensitivity.png", p, width = 4, height = 9)
ggsave("figure/Simulation_ATE_Sensitivity.pdf", p, width = 4, height = 9)






###'######################################################################
###'
###' Plot #2. ATEs and QTEs 
###'
###' 
###'

### Prepare dataset for plotting

df_AQTE <- tau_summary_collect 
df_AQTE$estimate <- "estimate"



### Plot!  

p <- ggplot(data = df_AQTE, aes(x = estimate, group = tau, col = tau)) +
  
  # Point
  geom_point(aes(y = mean), size = 3) + 
  
  # Error bar
  geom_errorbar(aes(ymin = pct2.5, ymax = pct97.5), width = 0.2) + 
  
  
  # Value label (posterior mean & sd)
  geom_text(aes(y = mean, label = sprintf("%.3f", round(mean, 3)), 
                hjust = 0.5, vjust = -1.0)) + 
  geom_text(aes(y = mean, label = paste0("(", sprintf("%.3f", round(sd, 3)), ")"), 
                hjust = -0.5, vjust = -1.0)) + 
  
  # Faceting
  facet_grid(rho ~ tau) + 
  
  # Flip coordinates
  coord_flip() + 
  
  
  # Theme settings
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank()) + 
  
  # Labels
  labs(title = "Estimated Treatment Effects under Bayesian Models",  
       subtitle = "Posterior Means, Standard Deviations (in parentheses), and 95% Credible Intervals",  
       y = "",
       x = "Estimated Treatment Effects")


### Save the plot 

setwd("~/Stan/20180204_Bayesian Model-Based Inference for RCT")

ggsave("figure/Simulation_QTE_Sensitivity.png", p, width = 14, height = 9)
ggsave("figure/Simulation_QTE_Sensitivity.pdf", p, width = 14, height = 9)

