###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS2
# Last Update: Oct 1, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse) # data cleaning and manipulation
require(magrittr)  # syntax
require(xtable)    # regression table output
require(perm)      # permutation tests
require(boot)      # bootstrapping
require(ggplot2)   # plots

set.seed(22)
select = dplyr::select
setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS2')

###############################################################################