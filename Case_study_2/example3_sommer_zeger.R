#################################################################
###                                         
### TUTORIAL ON STAN FOR CAUSAL INFERENCE
###
### Example 3: Sommer-Zeger example from Imbens & Rubin (1997)
###
### Data overview:
### y: survival
### z: random assignment
### d: vitamin supplement
###


############################################################
###
### SETUP
###

setwd("~/Desktop/Stan Tutorial")

library(rstan)
library(shinyStan)
set_cppo(mode = "fast")




####################################################################################
###
### LOAD DATA: SOMMER-ZEGER EXAMPLE FROM IMBENS AND RUBIN (1997)
### 


y <- c( rep(0, 74),
        rep(1, 11514),
        rep(0, 34),
        rep(1, 2385),
        rep(0, 12),
        rep(1, 9663) )

z <- c( rep(0, 74 + 11514), 
        rep(1, 34 + 2385 + 12 + 9663))

d <- c( rep(0, 74 + 11514 + 34 + 2385), 
        rep(1, 12 + 9663))



stan_data <- list(  N = length(y),
                    y = y,
                    z = z,
                    d = d)



####################################################################################
###
### RUN STAN MODEL WITH EXCLUSION RESTRICTION
###


stan_fit_ER <- stan( file = "example3_sommer_zeger_wER.stan", 
                     data = stan_data, 
                     iter = 1000, chains = 2)



####################################################################################
###
### ESTIMATE CACE WITH EXCLUSION RESTRICTION
###

## plot output
plot(stan_fit_ER)

## extract parameters 
eta_c1 <- extract(stan_fit_ER, pars = "eta_c1")$eta_c1
eta_c0 <- extract(stan_fit_ER, pars = "eta_c0")$eta_c0


## calculate treatment effect (in per-10,000 units)
cace <- (eta_c1 - eta_c0)*10^3

## plot CACE
hist(cace, main = "ITT for Compliers with Exclusion Restriction", 
     xlab = "ITT for Compliers (per 10,000)", col = "grey", border = "white",
     xlim = c(-4, 10), breaks = 30)

## plot median CACE
abline(v = median(cace), col = 2)

## plot 90% CI
quantile(cace, probs = c(0.05, 0.95))
abline(v = quantile(cace, probs = c(0.05, 0.95)), lty = 2, col = 2)






####################################################################################
###
### RUN STAN MODEL WITHOUT EXCLUSION RESTRICTION FOR NEVER TAKERS
###


stan_fit_noER <- stan( file = "example3_sommer_zeger_noER.stan", 
                       data = stan_data, 
                       iter = 1000, chains = 2)



####################################################################################
###
### ESTIMATE CACE AND NACE (WITHOUT EXCLUSION RESTRICTION)
###


## plot output
plot(stan_fit_noER)

## extract parameters 
eta_c1_no_er <- extract(stan_fit_noER, pars = "eta_c1")$eta_c1
eta_c0_no_er <- extract(stan_fit_noER, pars = "eta_c0")$eta_c0

eta_nt1 <- extract(stan_fit_noER, pars = "eta_nt1")$eta_nt1
eta_nt0 <- extract(stan_fit_noER, pars = "eta_nt0")$eta_nt0


## calculate treatment effect (in per-10,000 units)
cace_no_er <- (eta_c1_no_er - eta_c0_no_er)*10^3
nace <- (eta_nt1 - eta_nt0)*10^3


## plot CACE
hist(cace_no_er, main = "ITT for Compliers without Exclusion Restriction", 
     xlab = "ITT for Compliers (per 10,000)", col = "grey", border = "white",
     xlim = c(-4, 10), breaks = 30)

## plot median CACE
abline(v = median(cace_no_er), col = 2)

## plot 90% CI
quantile(cace_no_er, probs = c(0.05, 0.95))
abline(v = quantile(cace_no_er, probs = c(0.05, 0.95)), lty = 2, col = 2)



## plot NACE
hist(nace, main = "ITT for Never Takers", 
     xlab = "ITT for Never Takers (per 10,000)", col = "grey", border = "white",
     xlim = c(-25, 40), breaks = 30)

## plot median CACE
abline(v = median(nace), col = 2)

## plot 90% CI
quantile(nace, probs = c(0.05, 0.95))
abline(v = quantile(nace, probs = c(0.05, 0.95)), lty = 2, col = 2)


## plot joint distribution
plot(cace_no_er, nace, type = 'n',
     xlab = "ITT for Compliers", ylab = "ITT for Never Takers",
     xlim = c(-4, 10), ylim = c(-30, 30), bty = 'n')
abline(h = 0)
abline(v = 0)
points(cace_no_er, nace, pch = 20)



##################################################################
###
### COMPARE WITH AND WITHOUT EXCLUSION RESTRICTION
### 

boxplot(cace, cace_no_er,
        names = c("With Excl. Rest.", "Without Excl. Rest."),
        col = "grey",
        ylab = "ITT for Compliers")
