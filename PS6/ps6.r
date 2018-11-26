###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS6
# Last Update: Nov 26, 2018
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
require(rdrobust)  # RD stuff
require(rdlocrand) # RD stuff
require(rddensity) # RD stuff

options(scipen = 999)
setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS6')

# read in dataset
df <- read_csv('HeadStart.csv')

###############################################################################
# Question 2) The Effect of Head Start on Child Mortality
###############################################################################
# Q2.1 (a): RD Plots

# four plots: {evenly-spaced, quantile-spaced} X (IMSE-optimal,data var}
png('r/1_1a.png')
rdplot(df$mort_related_pre, df$povrate60, binselect = 'es',
       title = 'Evenly-spaced binning, IMSE-optimal',
       x.label = 'Poverty Rate (Normalized)',
       y.label = 'Average Mortality Rate (Pre)')
dev.off()

png('r/1_1b.png')
rdplot(df$mort_related_pre, df$povrate60, binselect = 'esmv',
       title = 'Evenly-spaced binning, variability-mimicking',
       x.label = 'Poverty Rate (Normalized)',
       y.label = 'Average Mortality Rate (Pre)')
dev.off()

png('r/1_1c.png')
rdplot(df$mort_related_pre, df$povrate60, binselect = 'qs',
       title = 'Quantile-spaced binning, IMSE-optimal',
       x.label = 'Poverty Rate (Normalized)',
       y.label = 'Average Mortality Rate (Pre)')
dev.off()

png('r/1_1d.png')
rdplot(df$mort_related_pre, df$povrate60, binselect = 'qsmv',
       title = 'Quantile-spaced binning, variability mimicking',
       x.label = 'Poverty Rate (Normalized)',
       y.label = 'Average Mortality Rate (Pre)')
dev.off()

###############################################################################
# Q2.1 (b): Falsification Tests
# falsification of RD design using three methods

# (i) histogram plots
df %<>% mutate(treated = povrate60 > 0)

p1_2i <- ggplot(df, aes(x = povrate60, fill = treated, color = treated)) +
  geom_histogram() + theme_minimal() +
  labs(title = 'Histogram by Treatment Type',
       y = 'Frequency', x = 'Poverty Rate (Normalized)')
ggsave('r/1_2i.png')

# (ii) binomial tests
rdwinselect(df$povrate60, wmin = .05, wstep = .05, nwindows = 100)

# (iii) continuity-in-density tests
rddensity(df$povrate60, all = T) %>% summary
rddensity(df$povrate60, all = T, h = c(5,5)) %>% summary
rddensity(df$povrate60, all = T, h = c(1,1)) %>% summary

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