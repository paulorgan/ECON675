###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS3
# Last Update: Oct 17, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse) # data cleaning and manipulation
require(magrittr)  # syntax
require(ggplot2)   # plots
require(xtable)    # tables for LaTeX
require(boot)      # bootstrapping
require(gmm)       # GMM estimation

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

# add variable to indicate having outcome data
df %<>% mutate(s = 1-dmissing)

# add log variable for us in regression [log(S_incomepc + 1)]
df %<>% mutate(log_SincpcP1 = log(S_incomepc + 1))

###############################################################################
# Q2.2b - feasible estimator

# since we are using the logistic CDF as our F, we have simplified instruments
# here we construct the moment condition 0=E[g*m] which is just
# (instruments)*(y_i - F(t\theta+x\gamma))
# function of theta and gamma (parameters) and t and X (data)
g_mcar <- function(beta, data){
  theta <- beta[1]
  gamma <- beta[2:4]
  # constructing the F(t\theta+x\gamma) term
  F_component <- plogis(data$dpisofirme * theta +
                          data$S_age * gamma[1] +
                          data$S_HHpeople * gamma[2] +
                          data$log_SincpcP1 * gamma[3])
  
  vars <- c('dpisofirme', 'S_age', 'S_HHpeople', 'log_SincpcP1')
  Zvars <- paste0('Z_', vars)
  
  for(i in 1:length(vars)){
    data[,Zvars[i]] <- data[,vars[i]]*(data$danemia - F_component)
  }
  
  out <- data[,Zvars] %>% as.matrix
  
  return(out)
}

# bootstrap statistic
boot_stat <- function(dat, i){
  gmm(g_mcar, dat[i,], t0=c(0,0,0,0), wmatrix='ident', vcov='iid')$coef
}

# some incomplete values seem to be messing things up (drops three observations)
df <- df[complete.cases(df[,c('dpisofirme','S_age','S_HHpeople','log_SincpcP1')]),]

# run bootstrap, 999 replications
ptm <- proc.time()
set.seed(22)
boot_results <- boot(data=df[df$s==1, ], R=999, statistic = boot_stat, stype = "i")
proc.time() - ptm
# runtime is 24 minutes

# generate table of desired results (showing same results as in Q1)
t3 <- matrix(NA, nrow=4, ncol=6) %>% as.data.frame
names(t3) <- c('Coefficient', 'StdError', 'tValue', 'pValue', 'CIlow', 'CIhigh')

t3$Coefficient <- boot_results$t0
t3$StdError    <- apply(boot_results$t, 2, sd)
t3$tValue      <- t3$Coefficient / t3$StdError

# define function to calculate pvalue from bootstrap results
getPval <- function(coef, coef_b){
  2 * max( mean( (coef_b-coef) >= abs(coef) ), mean( (coef_b-coef) <= -1*abs(coef) ) )
}

# apply over covariates, also get confidence intervals
for(i in 1:nrow(t3)){
  t3$pValue[i] <- getPval(boot_results$t0[i], boot_results$t[,i])
  
  t3$CIlow[i]  <- quantile(boot_results$t[,i], .025)
  t3$CIhigh[i] <- quantile(boot_results$t[,i], .975)
}

# add rownames
row.names(t3) <- c('dpisofirme', 'S_age', 'S_HHpeople', 'log(S_incomepc+1)')

# output for LaTeX
xtable(t3, digits = c(0,3,3,3,4,3,3))

###############################################################################
# Q2.3c - feasible estimator

# similar to above, but first estimate propensity score then plug into gmm
# only difference in this function is to add (*ipw) (weighting observations)
g_mar <- function(beta, data){
  # restrict to nonmissing data
  data <- data[data$s==1,]
  
  # beta has two components
  theta <- beta[1]
  gamma <- beta[2:4]
  
  # constructing the F(t\theta+x\gamma) term
  F_component <- plogis(data$dpisofirme * theta +
                          data$S_age * gamma[1] +
                          data$S_HHpeople * gamma[2] +
                          data$log_SincpcP1 * gamma[3])
  
  vars <- c('dpisofirme', 'S_age', 'S_HHpeople', 'log_SincpcP1')
  Zvars <- paste0('Z_', vars)
  
  for(i in 1:length(vars)){
    data[,Zvars[i]] <- data[,vars[i]]*(data$danemia - F_component)*data$ipw
  }
  
  out <- data[,Zvars] %>% as.matrix
  
  return(out)
}

# bootstrap statistic
# this has more changes - need to estimate p using logit
# then plug that into our data as "ipw" and run gmm
boot_stat <- function(dat, i){
  boot_data <- dat[i,]
  
  # estimate p using logit
  reg <- glm(s ~ dpisofirme + S_age + S_HHpeople + log_SincpcP1 - 1,
             data = boot_data, family = binomial(link = 'logit'))
  
  # get fitted values
  predicted_p <- reg$fitted
  
  # construct inverse prob weights
  boot_data$ipw <- 1 / predicted_p
  
  # run GMM using the constructed weights
  gmm(g_mar, boot_data, t0=c(0,0,0,0), wmatrix='ident', vcov='iid')$coef
}

# some incomplete values seem to be messing things up (drops three observations)
df <- df[complete.cases(df[,c('dpisofirme','S_age','S_HHpeople','log_SincpcP1')]),]

# run bootstrap, 999 replications
ptm <- proc.time()
set.seed(22)
boot_results <- boot(data=df, R=999, statistic = boot_stat, stype = "i")
proc.time() - ptm
# runtime is 35mins

# generate table of desired results (showing same results as in Q1)
t4 <- matrix(NA, nrow=4, ncol=6) %>% as.data.frame
names(t4) <- c('Coefficient', 'StdError', 'tValue', 'pValue', 'CIlow', 'CIhigh')

t4$Coefficient <- boot_results$t0
t4$StdError    <- apply(boot_results$t, 2, sd)
t4$tValue      <- t4$Coefficient / t4$StdError

# define function to calculate pvalue from bootstrap results
getPval <- function(coef, coef_b){
  2 * max( mean( (coef_b-coef) >= abs(coef) ), mean( (coef_b-coef) <= -1*abs(coef) ) )
}

# apply over covariates, also get confidence intervals
for(i in 1:nrow(t4)){
  t4$pValue[i] <- getPval(boot_results$t0[i], boot_results$t[,i])
  
  t4$CIlow[i]  <- quantile(boot_results$t[,i], .025)
  t4$CIhigh[i] <- quantile(boot_results$t[,i], .975)
}

# add rownames
row.names(t4) <- c('dpisofirme', 'S_age', 'S_HHpeople', 'log(S_incomepc+1)')

# output for LaTeX
xtable(t4, digits = c(0,3,3,3,4,3,3))

###############################################################################
# Q2.3d - feasible estimator, with trimming

# same as above, but trimming observations with p<.1 (1/p>10)
g_mar_trim <- function(beta, data){
  # restrict to nonmissing data
  data <- data[data$s==1,]
  
  # trim observations
  data <- data[data$ipw<=10,]
  
  # beta has two components
  theta <- beta[1]
  gamma <- beta[2:4]
  
  # constructing the F(t\theta+x\gamma) term
  F_component <- plogis(data$dpisofirme * theta +
                          data$S_age * gamma[1] +
                          data$S_HHpeople * gamma[2] +
                          data$log_SincpcP1 * gamma[3])
  
  vars <- c('dpisofirme', 'S_age', 'S_HHpeople', 'log_SincpcP1')
  Zvars <- paste0('Z_', vars)
  
  for(i in 1:length(vars)){
    data[,Zvars[i]] <- data[,vars[i]]*(data$danemia - F_component)*data$ipw
  }
  
  out <- data[,Zvars] %>% as.matrix
  
  return(out)
}

# bootstrap statistic
# only change is to use the g_mar_trim function instead of g_mar
boot_stat_trim <- function(dat, i){
  boot_data <- dat[i,]
  
  # estimate p using logit
  reg <- glm(s ~ dpisofirme + S_age + S_HHpeople + log_SincpcP1 - 1,
             data = boot_data, family = binomial(link = 'logit'))
  
  # get fitted values
  predicted_p <- reg$fitted
  
  # construct inverse prob weights
  boot_data$ipw <- 1 / predicted_p
  
  # run GMM using the constructed weights
  gmm(g_mar_trim, boot_data, t0=c(0,0,0,0), wmatrix='ident', vcov='iid')$coef
}

# run bootstrap, 999 replications
ptm <- proc.time()
set.seed(22)
boot_results <- boot(data=df, R=999, statistic = boot_stat_trim, stype = "i")
proc.time() - ptm
# runtime is 39mins

# generate table of desired results (showing same results as in Q1)
t5 <- matrix(NA, nrow=4, ncol=6) %>% as.data.frame
names(t5) <- c('Coefficient', 'StdError', 'tValue', 'pValue', 'CIlow', 'CIhigh')

t5$Coefficient <- boot_results$t0
t5$StdError    <- apply(boot_results$t, 2, sd)
t5$tValue      <- t5$Coefficient / t5$StdError

# apply over covariates, also get confidence intervals
for(i in 1:nrow(t5)){
  t5$pValue[i] <- getPval(boot_results$t0[i], boot_results$t[,i])
  
  t5$CIlow[i]  <- quantile(boot_results$t[,i], .025)
  t5$CIhigh[i] <- quantile(boot_results$t[,i], .975)
}

# add rownames
row.names(t5) <- c('dpisofirme', 'S_age', 'S_HHpeople', 'log(S_incomepc+1)')

# output for LaTeX
xtable(t5, digits = c(0,3,3,3,4,3,3))

###############################################################################
### Question 3: When Bootstrap Fails
###############################################################################
# clean up
rm(list = ls())
gc()

###############################################################################
# Q3.1 - nonparametric bootstrap

###############################################################################
# Q3.2 - parametric bootstrap

###############################################################################