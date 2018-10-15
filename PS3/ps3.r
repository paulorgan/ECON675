###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS3
# Last Update: Oct 15, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse) # data cleaning and manipulation
require(magrittr)  # syntax
require(ggplot2)   # plots
require(xtable)    # tables for LaTeX

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


# Q1.9b - nonparametric bootstrap

# Q1.9c - propensity scores

###############################################################################
### Question 2: Semiparametric GMM with Missing Data
###############################################################################
# Q2.2b - feasible estimator

# Q2.3c - feasible estimator

###############################################################################
### Question 3: When Bootstrap Fails
###############################################################################
# Q3.1 - nonparametric bootstrap

# Q3.2 - parametric bootstrap

###############################################################################