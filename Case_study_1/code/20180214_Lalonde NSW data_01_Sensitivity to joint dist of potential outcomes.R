
###'######################################################################
###'
###' Bayesian Model-Based Inference for Completely Randomized Experiments
###' 
###' Sensitivity analysis for 
###' the assumption of joint distribution of potential outcomes
###'
###'
###' 20180214 JoonHo Lee
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
###' Preparing for loop
###' 
###' 1) Loop over with/without covariates
###' 
###' 2) Loop over rho-s
###'
###'

### Define looping variables

covariates_vec <- c("without covariates", "with covariates")

rho_vec <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)



### Define an array for saving stanfits

stanfit_collect <- array(list(), dim = c(length(rho_vec), length(covariates_vec)))




###'######################################################################
###'
###' Implement for loops
###'
###'

for (i in seq(length(rho_vec))){  # (1) Loop over 6 rho-s
  
  for (j in seq(length(covariates_vec))) {  # (2) Loop over 2 covariates conditions
    
    
    ### Assign loop variables
    rho_value <- rho_vec[i]
    covariates_cond <- covariates_vec[j]
    
    
    ### Assign data list for stan
    stan_data <- list(  y = y,
                        z = z,
                        x = x,
                        y_pos = as.numeric(lalonde$re78 > 0), ## for two-part model
                        N = nrow(lalonde),
                        N_cov = ncol(x), 
                        rho = rho_value )
    
    ### Fit stan model
    
    if (rho_value == 1.0){  # perfect correlation (rho == 1.0)
      
      if (covariates_cond == "without covariates"){
        
        stan_fit <- stan(file = "code/Model_03_General Form without Covariates_Perfect correlation.stan",
                         data = stan_data,
                         iter = 1000, chains = 4)
      }else{  # "with covariates"
        stan_fit <- stan(file = "code/Model_04_General Form with Covariates_Perfect correlation.stan",
                         data = stan_data,
                         iter = 1000, chains = 4)
      }
    }else{  # non-perfect correlation (rho != 1.0)
      
      if (covariates_cond == "without covariates"){
        
        stan_fit <- stan(file = "code/Model_01_General Form without Covariates.stan",
                         data = stan_data,
                         iter = 1000, chains = 4)
      }else{  # "with covariates"
        stan_fit <- stan(file = "code/Model_02_General Form with Covariates.stan",
                         data = stan_data,
                         iter = 1000, chains = 4)
      }
    }
    
    
    ### Store stanfit in the array 
    stanfit_collect[[i, j]] <- stan_fit
    
    
  }  # End of loop over 2 covariates conditions
}    # End of loop over 6 rho-s



###'######################################################################
###'
###' Save the resulting array
###'
###'

save(stanfit_collect, file = "data/stanfit_collect.RData")




###'######################################################################
###'
###' Generate table
###' 
###' (1) Long data for plotting and tabulating
###'
###'


### Load stanfit collect

load(file = "data/stanfit_collect.RData")



###' Generate empty data frame to collect summary statistics 
###' for posterior distribution of tau

tau_summary_collect <- data.frame(matrix(ncol = 13, nrow = 0))

var_names <- c("rho", "covariates", "tau", 
               "mean", "se_mean", "sd", 
               "pct2.5", "pct25", "pct50", "pct75", "pct97.5", 
               "n_eff", "Rhat")

names(tau_summary_collect) <- var_names



### Start for loop to collect summary stats

for (i in seq(length(rho_vec))){  # (1) Loop over 6 rho-s
  
  for (j in seq(length(covariates_vec))) {  # (2) Loop over 2 covariates conditions
    
    
    ### Extract stanfit object
    
    stanfit <- stanfit_collect[[i, j]]
    
    
    ###' Assign summary statistics for the posterior distributions of
    ###' tau and tau_samp
    
    tau_summary <- summary(stanfit, c("tau", "tau_samp", 
                                      "tau_qte25", "tau_qte50", "tau_qte75"))$summary
    
    
    ### Add labels
    
    rho_lab <- c(paste0("rho = ", sprintf("%.1f", rho_vec[i])))
    
    covariates_lab <- covariates_vec[j]
    
    add_label <- matrix(rep(cbind(rho_lab, covariates_lab), 5), nrow = 5, byrow = TRUE)
    
    tau_df <- cbind.data.frame(add_label, 
                               rownames(tau_summary), 
                               tau_summary)
    
    
    ### Assign row names and column names
    rownames(tau_df) <- c()
    
    colnames(tau_df) <- var_names
    
    
    
    ### Append to the tau_summary_collect
    
    tau_summary_collect <- rbind.data.frame(tau_summary_collect, tau_df)
  }
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



### Add guidelines using estimates from Neyman's approach

summary(lm(y ~ z))  # Without covariates

summary(lm(y ~ x + z)) # With covariates



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
  facet_grid(rho ~ covariates) + 
  
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
  labs(title = "Estimated Treatment Effects under Bayesian Models for NSW Program Data",  
       subtitle = "Posterior Means, Standard Deviations (in parentheses), and 95% Credible Intervals",  
       y = "",
       x = "Estimated Average Treatment Effect")


### Save the plot 

setwd("~/Stan/20180204_Bayesian Model-Based Inference for RCT")

ggsave("figure/Lalonde_ATE_Sensitivity.png", p, width = 8, height = 9)
ggsave("figure/Lalonde_ATE_Sensitivity.pdf", p, width = 8, height = 9)






###'######################################################################
###'
###' Plot #2. ATEs and QTEs 
###'
###' 
###'

### Prepare dataset for plotting

df_AQTE <- tau_summary_collect 



### Plot!  

p <- ggplot(data = df_AQTE, aes(x = covariates, group = tau, col = covariates)) +
  
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
  
  
  # Color setting
  scale_colour_manual(labels = c("Without covariates", 
                                 "With covariates"), 
                      values = c("dodgerblue1", "firebrick1")) + 
  
  
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
  labs(title = "Estimated Treatment Effects under Bayesian Models for NSW Program Data",  
       subtitle = "Posterior Means, Standard Deviations (in parentheses), and 95% Credible Intervals",  
       y = "",
       x = "Estimated Treatment Effects")


### Save the plot 

setwd("~/Stan/20180204_Bayesian Model-Based Inference for RCT")

ggsave("figure/Lalonde_QTE_Sensitivity.png", p, width = 14, height = 9)
ggsave("figure/Lalonde_QTE_Sensitivity.pdf", p, width = 14, height = 9)

