###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS5
# Last Update: Oct 29, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse) # data cleaning and manipulation
require(magrittr)  # syntax
require(ggplot2)   # plots
require(sandwich)  # robust standard errors
require(ivpack)    # IV regression
require(xtable)    # tables for LaTeX
require(stargazer) # tables for LaTeX

options(scipen = 999)
setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS5')

###############################################################################
# Question 2) Weak Instrument Simulations
###############################################################################

# sample size
n <- 200

# set of gammas we will test
gammas <- c(0,0.25,9,99)
gammas <- gammas / n %>% sqrt()

# define data-generating process
dgp <- function(n, gamma){
  z <- rnorm(n, 0, 1)
  u <- rnorm(n, 0, 1)
  v <- .99*u + sqrt(1-.99^2)*rnorm(n, 0, 1)
  
  x <- gamma * z + v
  y <- u # beta is 0
  
  out = data.frame(y = y, x = x, z = z)
  return(out)
}

# testing for now
gamma <- gammas[2]

# define function to estimate everything we want
bigFunction <- function(rep, gamma){
  df <- dgp(n, gamma)
  
  # ols regression
  ols <- lm(y ~ x, data = df)
  
  # ols coef and se
  beta_ols <- ols$coefficients['x']
  se_ols   <- sqrt(vcovHC(ols,'HC1')['x','x'])
  
  # ols rejection indicator (null beta=0)
  rej_ols <- 1 * (beta_ols/se_ols > 1.96)
  
  # iv regression
  iv <- ivreg(y ~ x | z, data = df)
  
  # iv coef and se
  beta_iv <- iv$coefficients['x']
  se_iv   <- sqrt(vcovHC(iv,'HC1')['x','x'])
  
  # iv rejection indicator (null beta=0)
  rej_iv <- 1 * (beta_iv/se_iv > 1.96)
  
  # iv F stat
  F_iv <- summary(iv)$waldtest[1]
  
  # single row dataframe to report, which we can combine over replications
  out <- data.frame(rep = rep, gamma = gamma,
                    beta_ols=beta_ols, se_ols=se_ols, rej_ols=rej_ols,
                    beta_iv=beta_iv, se_iv=se_iv, rej_iv=rej_iv, F_iv=F_iv)
  return(out)
}

# run function M times for each gamma, save results
set.seed(22)

M    <- 5000
reps <- 1:M

ptm <- proc.time()
results1 <- lapply(reps, FUN = bigFunction, gamma=gammas[1]) %>% bind_rows
results2 <- lapply(reps, FUN = bigFunction, gamma=gammas[2]) %>% bind_rows
results3 <- lapply(reps, FUN = bigFunction, gamma=gammas[3]) %>% bind_rows
results4 <- lapply(reps, FUN = bigFunction, gamma=gammas[4]) %>% bind_rows
proc.time() - ptm
# runtime is 4 minutes

# summarize results of a single column
summarizeResult <- function(col){
  m  <- mean(col)
  sd <- sd(col)
  q1 <- quantile(col, 0.1)
  q5 <- quantile(col, 0.5)
  q9 <- quantile(col, 0.9)
  out = c(m, sd, q1, q5, q9)
  return(out)
}

# summarize results for a gamma simulation
summarizeResults <- function(data, gamma){
  df <- data
  
  # empty matrix to fill
  tab <- matrix(NA,7,5)
  
  tab[1,] <- summarizeResult(df$beta_ols)
  tab[2,] <- summarizeResult(df$se_ols)
  tab[3,] <- summarizeResult(df$rej_ols)
  tab[4,] <- summarizeResult(df$beta_iv)
  tab[5,] <- summarizeResult(df$se_iv)
  tab[6,] <- summarizeResult(df$rej_iv)
  tab[7,] <- summarizeResult(df$F_iv)

  # append gamma as a column to keep track
  out <- cbind(gamma, tab) %>% as.data.frame()
  names(out) <- c('gamma', 'mean', 'sd', 'q1', 'q5', 'q9')
  return(out)
}

# apply functions, combine, write to csv to make LaTeX table
summ1 <- summarizeResults(results1, gammas[1])
summ2 <- summarizeResults(results2, gammas[2])
summ3 <- summarizeResults(results3, gammas[3])
summ4 <- summarizeResults(results4, gammas[4])

summ <- bind_rows(summ1, summ2, summ3, summ4)
summ %>% write_csv('q2_r.csv')

###############################################################################
# Question 3: Weak Instrument - Empirical Study
###############################################################################

# clean up
rm(list = ls())
gc()

###############################################################################
# Question 3.1) Angrist and Krueger
df <- read_csv('Angrist_Krueger.csv')

# OLS1
ols1 <- lm(l_w_wage ~ educ + non_white + SMSA + married +
             factor(region) + factor(YoB_ld), data = df)

# 2SLS1
iv1 <- ivreg(l_w_wage ~ educ + non_white + SMSA + married +
               factor(region) + factor(YoB_ld) |
               factor(QoB)*factor(YoB_ld) + non_white + SMSA + married +
               factor(region) + factor(YoB_ld), data = df)

# OLS2
ols2 <- lm(l_w_wage ~ educ + non_white + SMSA + married + age_q + age_sq +
             factor(region) + factor(YoB_ld), data = df)

# 2SLS1
iv2 <- ivreg(l_w_wage ~ educ + non_white + SMSA + married + age_q + age_sq +
               factor(region) + factor(YoB_ld) |
               factor(QoB)*factor(YoB_ld) + non_white + SMSA + married +
               age_q + age_sq + factor(region) + factor(YoB_ld), data = df)

# output for LaTeX
vars <- c('educ', 'non_white', 'SMSA', 'married', 'age_q', 'age_sq')

stargazer(ols1, iv1, ols2, iv2, keep = vars, # show all four regs
          dep.var.labels.include = F, # drop the l_w_wage heading
          model.names = F, # we replace them below
          column.labels = c('OLS 1', '2SLS 1', 'OLS 2', '2SLS 2'),
          add.lines = list(c('9 Year-of-birth dummies', rep('Yes', 4)),
                           c('8 Region-of-residence dummies', rep('Yes', 4))),
          omit.table.layout = 'sn', # get rid of end table stuff (N, R2, etc.)
          star.char = c('', '', ''), # remove asterisks
          digits = 4) # decimals

###############################################################################
# Question 3.2) Bound, Jaeger, and Baker

###############################################################################