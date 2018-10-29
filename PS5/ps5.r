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
require(xtable)    # tables for LaTeX
require(stargazer) # tables for LaTeX
require(ivpack)    # IV regression

setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS5')

###############################################################################
# Question 2) Weak Instrument Simulations
###############################################################################


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