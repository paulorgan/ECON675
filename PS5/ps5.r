###############################################################################
# Author:  Paul Organ
# Date:    Sept 11, 2018
# Purpose: ECON 675, PS5
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

require(tidyverse)
require(magrittr)

# for table output
require(xtable)

select = dplyr::select

setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS5')

###############################################################################
# Question 2) Weak Instrument Simulations

###############################################################################
# Question 3.1) Angrist and Krueger
df <- read_csv('Angrist_Krueger.csv')

# OLS1
ak_1 <- lm(l_w_wage ~ educ + non_white + SMSA + married +
             factor(region) + factor(YoB_ld), df)

# OLS2
ak_2 <- lm(l_w_wage ~ educ + non_white + SMSA + married + age_q + age_sq +
             factor(region) + factor(YoB_ld), df)

# 2SLS1

###############################################################################
# Question 3.2) Bound, Jaeger, and Baker

###############################################################################