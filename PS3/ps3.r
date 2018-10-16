###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS3
# Last Update: Oct 16, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse) # data cleaning and manipulation
require(magrittr)  # syntax
require(ggplot2)   # plots
require(xtable)    # tables for LaTeX
require(boot)      # bootstrapping

setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS3')

###############################################################################
### Question 1: Non-linear Least Squares
###############################################################################
# Q1.9 setup

# read in data
df <- read_csv('pisofirme.csv')

# add variable to indicate having outcome data
df %<>% mutate(s = 1-dmissing)

# add log variable for us in regression [log(S_incomepc + 1)]
df %<>% mutate(log_SincpcP1 = log(S_incomepc + 1))

# logistic regression model
reg <- glm(s ~ S_age + S_HHpeople + log_SincpcP1,
           data = df, family = binomial(link = 'logit'))

# Q1.9a - estimates and statistics
# extract coef, var, t-stat, p-val, 95% CI
t1 <- xtable(reg)
names(t1) <- c('Coefficient', 'StdError', 'zValue', 'pValue')

# add confidence interval
t1$CIlow  <- t1$Coefficient - 1.96*t1$StdError
t1$CIhigh <- t1$Coefficient + 1.96*t1$StdError

# output to LaTeX
xtable(t1, digits = c(0,3,3,3,4,3,3))

###############################################################################
# Q1.9b - nonparametric bootstrap

# define statistic to bootstrap
boot_stat <- function(df, i){
  reg <- glm(s ~ S_age + S_HHpeople + log_SincpcP1,
             data = df[i, ], family = binomial(link = 'logit'))
  return(reg$coefficients)
}

# run bootstrap, 999 replications
set.seed(22)
boot_results <- boot(data = df, R = 999, statistic = boot_stat)

# generate table of desired results
t2 <- matrix(NA, nrow=4, ncol=6) %>% as.data.frame
names(t2) <- c('Coefficient', 'StdError', 'tValue', 'pValue', 'CIlow', 'CIhigh')

t2$Coefficient <- boot_results$t0
t2$StdError    <- apply(boot_results$t, 2, sd)
t2$tValue      <- t2$Coefficient / t2$StdError

# define function to calculate pvalue from bootstrap results
getPval <- function(coef, coef_b){
  2 * max( mean( (coef_b-coef) >= abs(coef) ), mean( (coef_b-coef) <= -1*abs(coef) ) )
}

# apply over covariates, also get confidence intervals
for(i in 1:nrow(t2)){
  t2$pValue[i] <- getPval(boot_results$t0[i], boot_results$t[,i])
  
  t2$CIlow[i]  <- quantile(boot_results$t[,i], .025)
  t2$CIhigh[i] <- quantile(boot_results$t[,i], .975)
}

# add rownames
row.names(t2) <- row.names(t1)

# output for LaTeX
xtable(t2, digits = c(0,3,3,3,4,3,3))

###############################################################################
# Q1.9c - propensity scores

# use original logistic regression
reg <- glm(s ~ S_age + S_HHpeople + log_SincpcP1,
           data = df, family = binomial(link = 'logit'))

# get data for use in plots
props <- reg$fitted.values %>% as.data.frame
names(props) <- 'fitted'

# plot with a variety of kernels
p1 <- ggplot(props, aes(x=fitted)) +
  geom_histogram(aes(y=..density..)) +
  geom_density(kernel = 'gaussian', size=1, aes(color='Gaussian')) +
  geom_density(kernel = 'epanechnikov', size=1, aes(color='Epanechnikov')) +
  geom_density(kernel = 'triangular', size=1, aes(color='Triangular')) +
  scale_color_manual(name = 'Kernel Densities',
                     values = c(Gaussian='red',
                                Epanechnikov='blue',
                                Triangular='green')) +
  theme_minimal() +
  labs(x='Fitted Probability of Observation', y='Density')
p1

ggsave('q1_9c_R.png')

###############################################################################
### Question 2: Semiparametric GMM with Missing Data
###############################################################################
# clean up
rm(list = ls())
gc()

# read in data
df <- read_csv('pisofirme.csv')

###############################################################################
# Q2.2b - feasible estimator


###############################################################################
# Q2.3c - feasible estimator

###############################################################################
### Question 3: When Bootstrap Fails
###############################################################################
# Q3.1 - nonparametric bootstrap

###############################################################################
# Q3.2 - parametric bootstrap

###############################################################################