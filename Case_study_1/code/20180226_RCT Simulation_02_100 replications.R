
###'######################################################################
###'
###' Bayesian Model-Based Inference for Completely Randomized Experiments
###' 
###' Simulating RCT data - With covariates
###' 
###' Simulating over 100 replications
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
###' Loop over 100 replicated datasets
###'

### Define looping variables

n_rep <- 100



### Define an array for saving stanfits

stanfit_collect <- array(list(), dim = n_rep)





###'######################################################################
###'
###' Simulate RCT Data With Covariates
###' 
###' Loop over 100 replicated datasets
###'
###'


for (i in seq(n_rep)){  # (1) Loop over 100 replicated datasets
  
  
  ###'#######################
  ###' Simulate RCT Data
  ###' Assumming rho = 0.0
  ###'#######################
  
  N <- 500          # number of observations
  
  betay <- c(0.5, 1.0, -1.5)     # coefficients of X in the outcome model
  
  tau <- 0.25     # treatment effect
  
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
                    rho = 0.0)
  
  
  
  ###'#######################
  ###' Fit the Stan model
  ###'#######################
  
  mod <- stan_model("code/Model_02_General Form with Covariates.stan")
  
  stan_fit <- sampling(mod,
                       data = stan_data,
                       iter = 1000, chains = 4)

  
  ### Store stanfit in the array 
  stanfit_collect[[i]] <- stan_fit
  
  
  ### Delete DLL files for stable looping
  dso_filenames <- dir(tempdir(), pattern =.Platform$dynlib.ext)
  filenames <- dir(tempdir())
  
  for (j in seq(dso_filenames))
    dyn.unload(file.path(tempdir(), dso_filenames[j]))
  
  
  ### Some files w/ long filenames that didn't like to be removeed
  for (j in seq(filenames))
    if (file.exists(file.path(tempdir(), filenames[j])) & nchar(filenames[j]) < 42) 
      file.remove(file.path(tempdir(), filenames[j]))
  
}    # End of loop over 100 datasets




###'######################################################################
###'
###' Save the resulting array
###'
###'

save(stanfit_collect, file = "data/stanfit_collect_simulation_rep100.RData")




###'######################################################################
###'
###' Generate table
###' 
###' Long data for plotting and tabulating
###'
###'


### Load stanfit collect

load(file = "data/stanfit_collect_simulation_rep100.RData")



###' Generate empty data frame to collect summary statistics 
###' for posterior distribution of tau

tau_summary_collect <- data.frame(matrix(ncol = 13, nrow = 0))

var_names <- c("iter", "tau", 
               "mean", "se_mean", "sd", 
               "pct2.5", "pct25", "pct50", "pct75", "pct97.5", 
               "n_eff", "Rhat")

names(tau_summary_collect) <- var_names



### Start for loop to collect summary stats

for (i in seq(n_rep)){  # (1) Loop over 100 datasets
  
  
  ### Extract stanfit object
  
  stanfit <- stanfit_collect[[i]]
  
  
  ###' Assign summary statistics for the posterior distributions of
  ###' tau and tau_samp
  
  tau_summary <- summary(stanfit, c("tau", "tau_samp", 
                                    "tau_qte25", "tau_qte50", "tau_qte75"))$summary
  
  
  ### Add labels
  
  nrep_lab <- i
  
  add_label <- matrix(rep(cbind(nrep_lab), 5), nrow = 5, byrow = TRUE)
  
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
###' Calculate performance evaluators and their Monte Carlo errors:
###' 
###' (1) Bias
###' (2) Relative error
###' (3) Coverage probability
###' 
###' Over 2 treatment effect estimates:
###'
###' (A) ATE finite sample
###' (B) ATE super population
###'
###'


### Prepare loop over two ATE estimates

tau_vec <- c("ATE Sample", "ATE Population")

var_names <- c("bias_est", "bias_mce", "MCSD_est", "MCSD_mce", 
               "SEbar_est", "SEbar_mce", "relerror_est", "relerror_mce", 
               "CovP_est", "CovP_mce")

evaluator_collect <- data.frame(matrix(ncol = length(tau_vec) + 1, nrow = length(var_names)))

colnames(evaluator_collect) <- c("estimates", tau_vec)

evaluator_collect$estimates <- var_names




### Start for loop 

for (i in seq(length(tau_vec))){
  
  
  ### Subset tau estimate collection
  
  df <- tau_summary_collect %>%
    filter(tau == tau_vec[i])
  
  
  ###'###################
  ###' Bias and MCSD
  ###'###################
  
  ### Point estimates
  
  bias <- mean(df$mean, na.rm = TRUE) - 0.25
  MCSD <- sd(df$mean, na.rm = TRUE)
  
  
  ### Monte Carlo errors
  
  bias_mce <- sqrt(MCSD^2/n_rep)
  MCSD_mce <- sqrt(MCSD^2/(2*(n_rep - 1)))
  
  
  
  ###'###################
  ###' Relative error
  ###'###################
  
  ### Model-based Standard Error (SE) = posterior standard deviation
  
  SE <- df$sd
  
  
  ### Average model-based SE
  
  SE2_bar <- mean(SE^2, na.rm = TRUE)
  SE_bar <- sqrt(SE2_bar)
  
  
  ### Variance of model-based SE
  
  varSE2 <- var(SE^2, na.rm = TRUE)
  
  
  ### Monte Carlo error of Average model-based SE
  
  # SE_bar_mce <- sqrt(varSE2/4*n_rep*SE2_bar) # Returns near zero values
  SE_bar_mce <- sd(SE, na.rm = TRUE)  # mean method

  
  ### Relative error and its Monte Carlo error
  
  relerror <- (SE_bar/MCSD) - 1
  relerror_mce <- (SE_bar/MCSD)*sqrt((varSE2/(4*n_rep*(SE_bar^4))) + (1/(2*(n_rep -1))))
  
  
  
  ###'###################
  ###' Coverage
  ###'###################
  
  ### 95% coverage probability
  
  CovP <- mean(1*(df$pct2.5 <= 0.25 & df$pct97.5 >= 0.25), na.rm = TRUE)
  
  
  ### Monte Carlo error of coverage probability
  
  CovP_mce <- sqrt(CovP*(1 - CovP)/n_rep)
  
  
  ### Collect performance evaluator estimates
  
  est_vec <- c(bias, bias_mce, MCSD, MCSD_mce, 
               SE_bar, SE_bar_mce, relerror, relerror_mce, 
               CovP, CovP_mce)
  
  
  ### Collect in predefined dataframe
  
  evaluator_collect[, i+1] <- round(est_vec, 7)
}

### Check out the resulting data frame

evaluator_collect




###'######################################################################
###'
###' Plot to compare ATE sample vs. ATE population
###'
###'

### Reshape data frame & separate a factor for the evaluators

df <- evaluator_collect %>%
  gather(tau, value, 2:3) %>% 
  separate(estimates, c("evaluator", "estimate")) %>%
  spread(estimate, value) 
  
df$evaluator <- factor(df$evaluator, 
                       levels = c("bias", "MCSD", "SEbar", "relerror", "CovP"), 
                       labels = c("Bias", 
                                  "MCSD", 
                                  "Average SE", 
                                  "Relative Error", 
                                  "Coverage Probability"))
df$tau <- factor(df$tau)

df <- df %>%
  arrange(evaluator, tau)


### Data for different lines for facets

vline_data <- data.frame(evaluator = c("Bias", "Relative Error"), 
                         vline = c(0, 0))

vline_data$evaluator <- factor(vline_data$evaluator)



### Plot with free scales

p <- ggplot(data = df, aes(y = tau, col = tau, #group = tau,
                           x = est, xmin = est - 1.96*mce, xmax = est + 1.96*mce)) +
  
  # Point
  geom_point(size = 3) + 
  
  # Error bar
  geom_errorbarh(height = 0.2) + 
  
  
  # Value label (posterior mean & sd)
  geom_text(aes(x = est, label = sprintf("%.3f", round(est, 3)),
                hjust = 0.5, vjust = -1.0)) +
  geom_text(aes(x = est, label = paste0("(", sprintf("%.3f", round(mce, 3)), ")"),
                hjust = -0.5, vjust = -1.0)) +
  
  
  # Reference line at est = 0
  geom_vline(data = vline_data, aes(xintercept = vline), 
             size = 0.5, col = "black", linetype = "longdash") +
  
  # Faceting
  facet_grid(. ~ evaluator, scales = "free") + 
  

  # Color setting
  scale_colour_manual(labels = c("Super-population ATE", "Finite-sample ATE"), 
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
  labs(title = "Performance Evaluators of ATE Estimators",  
       subtitle = "Estimates, Monte Carlo Errors (in parentheses), and 95% Confidence Intervals",  
       y = "",
       x = "Estimates")


### Save the plot 

setwd("~/Stan/20180204_Bayesian Model-Based Inference for RCT")

ggsave("figure/Performance_Evaluators.png", p, width = 14, height = 4)
ggsave("figure/Performance_Evaluators.pdf", p, width = 14, height = 4)



