###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS4
# Last Update: Oct 30, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse) # data cleaning and manipulation
require(magrittr)  # syntax
require(ggplot2)   # plots
require(sandwich)  # robust standard errors
require(xtable)    # tables for LaTeX
require(boot)      # bootstrapping
require(MatchIt)   # matching

setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS4')

###############################################################################
### Question 2: Estimating Average Treatment Effects
###############################################################################
# Q2 setup

# load data
df <- read_csv('LaLonde_all.csv')

# add treatment variable in expanded dataset
df %<>% mutate(t = ifelse(treat == 1, 1, 0))

# add other variables we will use in the analysis
df %<>% mutate(logre74 = log(re74 + 1), logre75 = log(re75 + 1),
               age2 = age^2, educ2 = educ^2, age3 = age^3,
               bXu74 = black * u74, eXre74 = educ * logre74)

# empty s to fill with results (1 + 6*3)
ate <- matrix(NA, nrow=19,ncol=8) %>% as.data.frame

# add headings
names(ate) <- c('e_tau', 'e_se', 'e_CI_l', 'e_CI_h', 'p_tau', 'p_se', 'p_CI_l', 'p_CI_h')

# we need one for ATE, one for ATT
att <- ate

# define three model specifications
va <- c('age', 'educ', 'black', 'hisp', 'married',
        'nodegr', 'logre74', 'logre75')
vb <- c(va, 'age2', 'educ2', 'u74', 'u75')
vc <- c(vb, 'age3', 'bXu74', 'eXre74')

a <- paste('re78', paste(c('t', va),
                         collapse = ' + '), sep = ' ~ ') %>% as.formula()
b <- paste('re78', paste(c('t', vb),
                         collapse = ' + '), sep = ' ~ ') %>% as.formula()
c <- paste('re78', paste(c('t', vc),
                         collapse = ' + '), sep = ' ~ ') %>% as.formula()

ta <- paste('t', paste(va, collapse = ' + '), sep = ' ~ ') %>% as.formula()
tb <- paste('t', paste(vb, collapse = ' + '), sep = ' ~ ') %>% as.formula()
tc <- paste('t', paste(vc, collapse = ' + '), sep = ' ~ ') %>% as.formula()

ya <- paste('re78', paste(va, collapse = ' + '), sep = ' ~ ') %>% as.formula()
yb <- paste('re78', paste(vb, collapse = ' + '), sep = ' ~ ') %>% as.formula()
yc <- paste('re78', paste(vc, collapse = ' + '), sep = ' ~ ') %>% as.formula()

###############################################################################
# Q2.1: Difference-in-Means

# experimental data
r1e <- lm(re78 ~ t, data = df %>% filter(treat != 2))

# extract tau
ate$e_tau[1]  <- r1e$coefficients['t']

# hederosketastic robust standard errors
ate$e_se[1]   <- diag(vcovHC(r1e, type = "HC2")) %>% sqrt() %>% .['t']

# PSID data
r1p <- lm(re78 ~ t, data = df %>% filter(treat != 0))

# extract tau
ate$p_tau[1]  <- r1p$coefficients['t']

# hederosketastic robust standard errors
ate$p_se[1]   <- diag(vcovHC(r1p, type = "HC2")) %>% sqrt() %>% .['t']

# ATE and ATT are same here
att[1,] <- ate[1,]

###############################################################################
# Q2.2: Linear Least-Squares

# three specifications for each

# OLS for experimental data
r2ea <- lm(a, data = df %>% filter(treat != 2))
r2eb <- lm(b, data = df %>% filter(treat != 2))
r2ec <- lm(c, data = df %>% filter(treat != 2))

# extract tau
ate$e_tau[2]  <- r2ea$coefficients['t']
ate$e_tau[3]  <- r2eb$coefficients['t']
ate$e_tau[4]  <- r2ec$coefficients['t']

# robust standard errors 
ate$e_se[2]   <- diag(vcovHC(r2ea, type = "HC1")) %>% sqrt() %>% .['t']
ate$e_se[3]   <- diag(vcovHC(r2eb, type = "HC1")) %>% sqrt() %>% .['t']
ate$e_se[4]   <- diag(vcovHC(r2ec, type = "HC1")) %>% sqrt() %>% .['t']

# OLS for PSID data
r2pa <- lm(a, data = df %>% filter(treat != 0))
r2pb <- lm(b, data = df %>% filter(treat != 0))
r2pc <- lm(c, data = df %>% filter(treat != 0))

# extract tau
ate$p_tau[2]  <- r2pa$coefficients['t']
ate$p_tau[3]  <- r2pb$coefficients['t']
ate$p_tau[4]  <- r2pc$coefficients['t']

# robust standard errors 
ate$p_se[2]   <- diag(vcovHC(r2pa, type = "HC1")) %>% sqrt() %>% .['t']
ate$p_se[3]   <- diag(vcovHC(r2pb, type = "HC1")) %>% sqrt() %>% .['t']
ate$p_se[4]   <- diag(vcovHC(r2pc, type = "HC1")) %>% sqrt() %>% .['t']

# again, ATE is same as ATT here
att[2:4,] <- ate[2:4,]

###############################################################################
# Q2.3: Regression Imputation (rows 5-7)

df_e <- df %>% filter(treat != 2)
df_p <- df %>% filter(treat != 0)

n_e <- nrow(df_e)
n_p <- nrow(df_p)

## covariates a
# treated
beta1 <- lm(ya, data = df %>% filter(treat == 1)) %>% .$coefficients
# untreated experimental
beta0_e <- lm(ya, data = df %>% filter(treat == 0)) %>% .$coefficients
# untreated PSID
beta0_p <- lm(ya, data = df %>% filter(treat == 2)) %>% .$coefficients

# tauhat for each
Z_e <- as.matrix(df_e[,va]) %>% cbind(1,.)
tauhat_e <- Z_e %*% (beta1-beta0_e)

Z_p <- as.matrix(df_p[,va]) %>% cbind(1,.)
tauhat_p <- Z_p %*% (beta1-beta0_p)

# save results (ATE)
ate[5,1] <- mean(tauhat_e)
ate[5,2] <- sd(tauhat_e)/sqrt(n_e) # not right!

ate[5,5] <- mean(tauhat_p)
ate[5,6] <- sd(tauhat_p)/sqrt(n_p) # not right!

# save results (ATT)
att[5,1] <- mean(tauhat_e[df_e$t==1])
att[5,2] <- sd(tauhat_e[df_e$t==1])/sqrt(sum(df_e$t==1)) # not right!

att[5,5] <- mean(tauhat_p[df_p$t==1])
att[5,6] <- sd(tauhat_p[df_p$t==1])/sqrt(sum(df_p$t==1)) # not right!

## covariates b
# treated
beta1 <- lm(yb, data = df %>% filter(treat == 1)) %>% .$coefficients
# untreated experimental
beta0_e <- lm(yb, data = df %>% filter(treat == 0)) %>% .$coefficients
# untreated PSID
beta0_p <- lm(yb, data = df %>% filter(treat == 2)) %>% .$coefficients

# tauhat for each
Z_e <- as.matrix(df_e[,vb]) %>% cbind(1,.)
tauhat_e <- Z_e %*% (beta1-beta0_e)

Z_p <- as.matrix(df_p[,vb]) %>% cbind(1,.)
tauhat_p <- Z_p %*% (beta1-beta0_p)

# save results (ATE)
ate[6,1] <- mean(tauhat_e)
ate[6,2] <- sd(tauhat_e)/sqrt(n_e) # not right!

ate[6,5] <- mean(tauhat_p)
ate[6,6] <- sd(tauhat_p)/sqrt(n_p) # not right!

# save results (ATT)
att[6,1] <- mean(tauhat_e[df_e$t==1])
att[6,2] <- sd(tauhat_e[df_e$t==1])/sqrt(sum(df_e$t==1)) # not right!

att[6,5] <- mean(tauhat_p[df_p$t==1])
att[6,6] <- sd(tauhat_p[df_p$t==1])/sqrt(sum(df_p$t==1)) # not right!

## covariates c
# treated
beta1 <- lm(yc, data = df %>% filter(treat == 1)) %>% .$coefficients
# untreated experimental
beta0_e <- lm(yc, data = df %>% filter(treat == 0)) %>% .$coefficients
# untreated PSID
beta0_p <- lm(yc, data = df %>% filter(treat == 2)) %>% .$coefficients

# tauhat for each
Z_e <- as.matrix(df_e[,vc]) %>% cbind(1,.)
tauhat_e <- Z_e %*% (beta1-beta0_e)

Z_p <- as.matrix(df_p[,vc]) %>% cbind(1,.)
tauhat_p <- Z_p %*% (beta1-beta0_p)

# save results (ATE)
ate[7,1] <- mean(tauhat_e)
ate[7,2] <- sd(tauhat_e)/sqrt(n_e) # not right!

ate[7,5] <- mean(tauhat_p)
ate[7,6] <- sd(tauhat_p)/sqrt(n_p) # not right!

# save results (ATT)
att[7,1] <- mean(tauhat_e[df_e$t==1])
att[7,2] <- sd(tauhat_e[df_e$t==1])/sqrt(sum(df_e$t==1)) # not right!

att[7,5] <- mean(tauhat_p[df_p$t==1])
att[7,6] <- sd(tauhat_p[df_p$t==1])/sqrt(sum(df_p$t==1)) # not right!

###############################################################################
# Q2.4: Inverse Probability Weighting (rows 8-10)

# trimming for PSID data drops a lot of observations
# not yet matching my Stata results

# first estimate propensity scores (we use a probit here)
r4ae <- glm(ta, data = df_e, family = binomial(link = 'probit'))
r4be <- glm(tb, data = df_e, family = binomial(link = 'probit'))
r4ce <- glm(tc, data = df_e, family = binomial(link = 'probit'))
df_e %<>% mutate(prop_a = predict(r4ae, df_e, type='response'),
                 prop_b = predict(r4be, df_e, type='response'),
                 prop_c = predict(r4ce, df_e, type='response'))

r4ap <- glm(ta, data = df_p, family = binomial(link = 'probit'))
r4bp <- glm(tb, data = df_p, family = binomial(link = 'probit'))
r4cp <- glm(tc, data = df_p, family = binomial(link = 'probit'))
df_p %<>% mutate(prop_a = predict(r4ap, df_p, type='response'),
                 prop_b = predict(r4bp, df_p, type='response'),
                 prop_c = predict(r4cp, df_p, type='response'))

# trim df_p for interior predicted propensities
df_p %<>% filter(prop_a >= .001, prop_b >= .001, prop_c >= .001,
                 prop_a <= .999, prop_b <= .999, prop_c <= .999)

# calculate probability of treatment
n_p <- nrow(df_p)

phat_e = sum(df_e$t==1)/n_e
phat_p = sum(df_p$t==1)/n_p

# build components for ATE
df_e %<>% mutate(left_a = t*re78/prop_a, right_a = (1-t)*re78/(1-prop_a),
                 left_b = t*re78/prop_b, right_b = (1-t)*re78/(1-prop_b),
                 left_c = t*re78/prop_c, right_c = (1-t)*re78/(1-prop_c))
df_p %<>% mutate(left_a = t*re78/prop_a, right_a = (1-t)*re78/(1-prop_a),
                 left_b = t*re78/prop_b, right_b = (1-t)*re78/(1-prop_b),
                 left_c = t*re78/prop_c, right_c = (1-t)*re78/(1-prop_c))

# calculate ATEs
ate[8,1]  <- mean(df_e$left_a) - mean(df_e$right_a)
ate[8,5]  <- mean(df_p$left_a) - mean(df_p$right_a)
ate[9,1]  <- mean(df_e$left_b) - mean(df_e$right_b)
ate[9,5]  <- mean(df_p$left_b) - mean(df_p$right_b)
ate[10,1] <- mean(df_e$left_c) - mean(df_e$right_c)
ate[10,5] <- mean(df_p$left_c) - mean(df_p$right_c)

# build components for ATT
df_e %<>% mutate(left_a = t*re78/phat_e, right_a = (prop_a/phat_e)*(1-t)*re78/(1-prop_a),
                 left_b = t*re78/phat_e, right_b = (prop_b/phat_e)*(1-t)*re78/(1-prop_b),
                 left_c = t*re78/phat_e, right_c = (prop_c/phat_e)*(1-t)*re78/(1-prop_c))
df_p %<>% mutate(left_a = t*re78/phat_p, right_a = (prop_a/phat_p)*(1-t)*re78/(1-prop_a),
                 left_b = t*re78/phat_p, right_b = (prop_b/phat_p)*(1-t)*re78/(1-prop_b),
                 left_c = t*re78/phat_p, right_c = (prop_c/phat_p)*(1-t)*re78/(1-prop_c))

# calculate ATTs
att[8,1]  <- mean(df_e$left_a) - mean(df_e$right_a)
att[8,5]  <- mean(df_p$left_a) - mean(df_p$right_a)
att[9,1]  <- mean(df_e$left_b) - mean(df_e$right_b)
att[9,5]  <- mean(df_p$left_b) - mean(df_p$right_b)
att[10,1] <- mean(df_e$left_c) - mean(df_e$right_c)
att[10,5] <- mean(df_p$left_c) - mean(df_p$right_c)

###############################################################################
# Q2.5: Doubly Robust

###############################################################################
# Q2.6: Nearest Neighbor Matching (rows 14-16)

# useful sources:
# https://cran.r-project.org/web/packages/MatchIt/vignettes/matchit.pdf

# redefine dataframes
df_e <- df %>% filter(treat != 2)
df_p <- df %>% filter(treat != 0)

# see https://r.iq.harvard.edu/docs/matchit/2.4-20/Examples2.html
# function to calculate ate and att using NN methods
# can't figure out how to pass a dataframe to the matchit command
# tried a million ways but none working, so writing two functions...
nnResultsE <- function(tmodel, ymodel, est){
    # using MatchIt
  reg_match <- matchit(tmodel, df_e, method = 'nearest')
  
  df_t <- match.data(reg_match, 'treat')
  df_c <- match.data(reg_match, 'control')
  
  # predict values and treatment effects
  reg_t       <- lm(ymodel, df_t)
  df_c$y_pred <- predict(reg_t, df_c)
  
  df_c %<>% mutate(teffect = y_pred-re78)
  
  reg_c <- lm(ymodel, df_c)
  df_t$y_pred  <- predict(reg_c, df_t)
  
  df_t %<>% mutate(teffect = re78-y_pred)
  
  # stack and predict treatment effects, etc.
  comb <- bind_rows(df_c, df_t)
  
  ate <- mean(comb$teffect)
  att <- mean(comb$teffect[comb$t==1])
  
  n_att = nrow(df_t)
  n_ate = n_att*2
  
  ate_se <- sd(comb$teffect)/sqrt(n_ate)
  att_se <- sd(comb$teffect[comb$t==1])/sqrt(att)
  
  out_ate = matrix(c(ate, ate_se, NA, NA),1,4)
  out_att = matrix(c(att, att_se, NA, NA),1,4)
  
  if(est == 'ate'){
    return(out_ate)
  } else {
    return(out_att)
  }
}

nnResultsP <- function(tmodel, ymodel, est){
  # using MatchIt
  reg_match <- matchit(tmodel, df_p, method = 'nearest')
  
  df_t <- match.data(reg_match, 'treat')
  df_c <- match.data(reg_match, 'control')
  
  # predict values and treatment effects
  reg_t       <- lm(ymodel, df_t)
  df_c$y_pred <- predict(reg_t, df_c)
  
  df_c %<>% mutate(teffect = y_pred-re78)
  
  reg_c <- lm(ymodel, df_c)
  df_t$y_pred  <- predict(reg_c, df_t)
  
  df_t %<>% mutate(teffect = re78-y_pred)
  
  # stack and predict treatment effects, etc.
  comb <- bind_rows(df_c, df_t)
  
  ate <- mean(comb$teffect)
  att <- mean(comb$teffect[comb$t==1])
  
  n_att = nrow(df_t)
  n_ate = n_att*2
  
  ate_se <- sd(comb$teffect)/sqrt(n_ate)
  att_se <- sd(comb$teffect[comb$t==1])/sqrt(att)
  
  out_ate = matrix(c(ate, ate_se, NA, NA),1,4)
  out_att = matrix(c(att, att_se, NA, NA),1,4)
  
  if(est == 'ate'){
    return(out_ate)
  } else {
    return(out_att)
  }
}

# save results
ate[14,1:4] <- nnResultsE(ta, ya, 'ate')
ate[14,5:8] <- nnResultsP(ta, ya, 'ate')

ate[15,1:4] <- nnResultsE(tb, yb, 'ate')
ate[15,5:8] <- nnResultsP(tb, yb, 'ate')

ate[16,1:4] <- nnResultsE(tc, yc, 'ate')
ate[16,5:8] <- nnResultsP(tc, yc, 'ate')

att[14,1:4] <- nnResultsE(ta, ya, 'att')
att[14,5:8] <- nnResultsP(ta, ya, 'att')

att[15,1:4] <- nnResultsE(tb, yb, 'att')
att[15,5:8] <- nnResultsP(tb, yb, 'att')

att[16,1:4] <- nnResultsE(tc, yc, 'att')
att[16,5:8] <- nnResultsP(tc, yc, 'att')

###############################################################################
# Q2.7: Propensity Score Matching (rows 17-19)

# using same code as above, with probit for distance (propensity score)
# https://www.r-bloggers.com/using-the-r-matchit-package-for-propensity-score-analysis/

# redefine dataframes
df_e <- df %>% filter(treat != 2)
df_p <- df %>% filter(treat != 0)

# only difference here is the matchit command uses probit
psResultsE <- function(tmodel, ymodel, est){
  # using MatchIt
  reg_match <- matchit(tmodel, df_e, method = 'nearest', distance = 'probit')
  
  df_t <- match.data(reg_match, 'treat')
  df_c <- match.data(reg_match, 'control')
  
  # predict values and treatment effects
  reg_t       <- lm(ymodel, df_t)
  df_c$y_pred <- predict(reg_t, df_c)
  
  df_c %<>% mutate(teffect = y_pred-re78)
  
  reg_c <- lm(ymodel, df_c)
  df_t$y_pred  <- predict(reg_c, df_t)
  
  df_t %<>% mutate(teffect = re78-y_pred)
  
  # stack and predict treatment effects, etc.
  comb <- bind_rows(df_c, df_t)
  
  ate <- mean(comb$teffect)
  att <- mean(comb$teffect[comb$t==1])
  
  n_att = nrow(df_t)
  n_ate = n_att*2
  
  ate_se <- sd(comb$teffect)/sqrt(n_ate)
  att_se <- sd(comb$teffect[comb$t==1])/sqrt(att)
  
  out_ate = matrix(c(ate, ate_se, NA, NA),1,4)
  out_att = matrix(c(att, att_se, NA, NA),1,4)
  
  if(est == 'ate'){
    return(out_ate)
  } else {
    return(out_att)
  }
}

psResultsP <- function(tmodel, ymodel, est){
  # using MatchIt
  reg_match <- matchit(tmodel, df_p, method = 'nearest', distance = 'probit')
  
  df_t <- match.data(reg_match, 'treat')
  df_c <- match.data(reg_match, 'control')
  
  # predict values and treatment effects
  reg_t       <- lm(ymodel, df_t)
  df_c$y_pred <- predict(reg_t, df_c)
  
  df_c %<>% mutate(teffect = y_pred-re78)
  
  reg_c <- lm(ymodel, df_c)
  df_t$y_pred  <- predict(reg_c, df_t)
  
  df_t %<>% mutate(teffect = re78-y_pred)
  
  # stack and predict treatment effects, etc.
  comb <- bind_rows(df_c, df_t)
  
  ate <- mean(comb$teffect)
  att <- mean(comb$teffect[comb$t==1])
  
  n_att = nrow(df_t)
  n_ate = n_att*2
  
  ate_se <- sd(comb$teffect)/sqrt(n_ate)
  att_se <- sd(comb$teffect[comb$t==1])/sqrt(att)
  
  out_ate = matrix(c(ate, ate_se, NA, NA),1,4)
  out_att = matrix(c(att, att_se, NA, NA),1,4)
  
  if(est == 'ate'){
    return(out_ate)
  } else {
    return(out_att)
  }
}

# save results
ate[17,1:4] <- psResultsE(ta, ya, 'ate')
ate[17,5:8] <- psResultsP(ta, ya, 'ate')

ate[18,1:4] <- psResultsE(tb, yb, 'ate')
ate[18,5:8] <- psResultsP(tb, yb, 'ate')

ate[19,1:4] <- psResultsE(tc, yc, 'ate')
ate[19,5:8] <- psResultsP(tc, yc, 'ate')

att[17,1:4] <- psResultsE(ta, ya, 'att')
att[17,5:8] <- psResultsP(ta, ya, 'att')

att[18,1:4] <- psResultsE(tb, yb, 'att')
att[18,5:8] <- psResultsP(tb, yb, 'att')

att[19,1:4] <- psResultsE(tc, yc, 'att')
att[19,5:8] <- psResultsP(tc, yc, 'att')

###############################################################################
# Q2 final steps

# confidence intervals
ate %<>% mutate(e_CI_l = e_tau - 1.96*e_se,
                e_CI_h = e_tau + 1.96*e_se,
                p_CI_l = p_tau - 1.96*p_se,
                p_CI_h = p_tau + 1.96*p_se)

att %<>% mutate(e_CI_l = e_tau - 1.96*e_se,
                e_CI_h = e_tau + 1.96*e_se,
                p_CI_l = p_tau - 1.96*p_se,
                p_CI_h = p_tau + 1.96*p_se)

rownames(ate) <- c(
  'Diff-in-Means',
  'OLS: a', 'OLS: b', 'OLS: c',
  'RI: a', 'RI: b', 'RI: c',
  'IPW: a', 'IPW: b', 'IPW: c',
  'DR: a', 'DR: b', 'DR: c',
  'NN: a', 'NN: b', 'NN: c',
  'Propensity: a', 'Propensity: b', 'Propensity: c'
)

rownames(att) <- rownames(ate)

# write to LaTeX
xtable(ate)
xtable(att)

###############################################################################
### Question 3: Post-model Selection Inference
###############################################################################
# Q3 setup
rm(list = ls())
gc()

###############################################################################
# Q3.1: Summary Statistics and Kernel Density Plots

# replications
M <- 1000

# sample size
n <- 50

# define data generating process
dgp <- function(n){
  x   <- rnorm(n,0,1)
  z   <- .85*x + sqrt(1-.85)*rnorm(n,0,1)
  eps <- rnorm(n,0,1)
  
  y = 1 + .5*x + z + eps
  
  out = data.frame(y=y,x=x,z=z)
  
  return(out)
}

# empty matrices to store estimates and indicator of coverage
est <- matrix(NA, nrow=M, ncol = 3)
cov <- matrix(NA, nrow=M, ncol = 3)

# replicate M times and save estimates
set.seed(22)
for(i in 1:M){
  # generate sample
  df <- dgp(n)
  
  # run 'long' regression: Equation (2)
  long <- lm(y ~ x + z, data = df)
  
  # extract first estimate
  beta_hat <- long$coefficients['x']
  est[i,1] <- beta_hat
  
  # get se to calculate confidence interval and check coverage
  se_hat <- sqrt(vcovHC(long, 'HC1')['x','x'])
  lb_hat <- beta_hat - 1.96*se_hat
  ub_hat <- beta_hat + 1.96*se_hat
  cov[i,1] <- lb_hat <= .5 & ub_hat >= .5
  
  # save gamma over se gamma (using robust std errors)
  gamma_hat <- long$coefficients['z']
  gamma_se  <- sqrt(vcovHC(long, 'HC1')['z','z'])
  tstat <- gamma_hat/gamma_se
  
  # run 'short' regression: Equation (3)
  short <- lm(y ~ x, data = df)
  
  # extract second estimate
  beta_tilde <- short$coefficients['x']
  est[i,2]   <- beta_tilde
  
  # get se to calculate confidence interval
  se_tilde <- sqrt(vcovHC(short, 'HC1')['x','x'])
  lb_tilde <- beta_tilde - 1.96*se_hat
  ub_tilde <- beta_tilde + 1.96*se_hat
  cov[i,2] <- lb_tilde <= .5 & ub_tilde >= .5
  
  # save third estimate (beta_check) based on t-stat on gamma
  est[i,3] <- ifelse(tstat >= 1.96, beta_hat, beta_tilde)
  cov[i,3] <- ifelse(tstat >= 1.96, cov[i,1], cov[i,2])
}

# summary statistic of distribution of each of the estimators
tab <- matrix(NA,nrow=4,ncol=3)
for(i in 1:3){
  tab[1,i] <- min(est[,i])
  tab[2,i] <- mean(est[,i])
  tab[3,i] <- median(est[,i])
  tab[4,i] <- max(est[,i])
}

# add labels
rownames(tab) <- c('min', 'mean', 'median', 'max')
colnames(tab) <- c('beta_hat', 'beta_tilde', 'beta_check')

# write to LaTeX
xtable(tab %>% round(3))

# kernel density plots
est %<>% as.data.frame
names(est) <- c('beta_hat', 'beta_tilde', 'beta_check')
est <- gather(est, key = 'Estimator', value = 'x')

# plot kernel densities
p <- ggplot(est, aes(x = x, group = Estimator, color = Estimator)) +
  geom_density(kernel = 'epanechnikov', size = 1) +
  theme_minimal() +
  labs(x = '', y = 'Density')
p

# save for LaTeX
ggsave('q3_1_R.png')

###############################################################################
# Q3.2: Empirical Coverage Rates

# we know beta_0 is actually .5, how often do we cover that?
# in the loop above, we calculated coverage for each replication
# now, summarize for each estimator:

cov_hat   <- mean(cov[,1])
cov_tilde <- mean(cov[,2])
cov_check <- mean(cov[,3])

###############################################################################