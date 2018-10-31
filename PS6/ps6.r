###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS6
# Last Update: Oct 31, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse) # data cleaning and manipulation
require(magrittr)  # syntax
require(ggplot2)   # plots
require(sandwich)  # robust standard errors
require(xtable)    # tables for LaTeX
require(stargazer) # tables for LaTeX

options(scipen = 999)
setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS6')

###############################################################################
# Question 2) The Effect of Head Start on Child Mortality
###############################################################################
# Q2.1 (a): RD Plots

# four plots: {evenly-spaced, quantile-spaced} X (IMSE-optimal,data var}

###############################################################################
# Q2.1 (b): Falsification Tests
# falsification of RD design using three methods

# (1) histogram plots

# (2) binomial tests

# (3) continuity-in-density tests

###############################################################################
# Q2.2 (a): constant treatment effect model

###############################################################################
# Q2.2 (b): heterogeneous treatment effect model

###############################################################################
# Q2.2 (c): local parametric model

###############################################################################
# Q2.3 (a): MSE-optimal RD estimators

###############################################################################
# Q2.3 (b): Robustness checks

###############################################################################
# Q2.4: Local Randomization Methods

###############################################################################