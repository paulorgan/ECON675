###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS6
# Last Update: Nov 27, 2018
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

source('multiplot.R') # multi-panel plots

options(scipen = 999)
setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS6')

# read in dataset
df <- read_csv('HeadStart.csv')

# rename for ease of use
df %<>% rename(pov = povrate60,
               rel_post = mort_related_post,
               rel_pre  = mort_related_pre,
               inj_post = mort_injury_post)

###############################################################################
# Question 2) The Effect of Head Start on Child Mortality
###############################################################################
# Q2.1.1: RD Plots

# four plots: {evenly-spaced, quantile-spaced} X (IMSE-optimal,data var}
png('r/1_1a.png')
rdplot(df$rel_pre, df$pov, binselect = 'es',
       title = 'Evenly-spaced binning, IMSE-optimal',
       x.label = 'Poverty Rate (Normalized)',
       y.label = 'Average Mortality Rate (Pre)')
dev.off()

png('r/1_1b.png')
rdplot(df$rel_pre, df$pov, binselect = 'esmv',
       title = 'Evenly-spaced binning, variability-mimicking',
       x.label = 'Poverty Rate (Normalized)',
       y.label = 'Average Mortality Rate (Pre)')
dev.off()

png('r/1_1c.png')
rdplot(df$rel_pre, df$pov, binselect = 'qs',
       title = 'Quantile-spaced binning, IMSE-optimal',
       x.label = 'Poverty Rate (Normalized)',
       y.label = 'Average Mortality Rate (Pre)')
dev.off()

png('r/1_1d.png')
rdplot(df$rel_pre, df$pov, binselect = 'qsmv',
       title = 'Quantile-spaced binning, variability mimicking',
       x.label = 'Poverty Rate (Normalized)',
       y.label = 'Average Mortality Rate (Pre)')
dev.off()

###############################################################################
# Q2.1.2: Falsification Tests
# falsification of RD design using three methods

# (i) histogram plots
df %<>% mutate(t = (pov > 0) * 1)

p1_2i <- ggplot(df, aes(x = pov, fill = t, color = treated)) +
  geom_histogram() + theme_minimal() +
  labs(title = 'Histogram by Treatment Type',
       y = 'Frequency', x = 'Poverty Rate (Normalized)')
ggsave('r/1_2i.png')

# (ii) binomial tests
rdwinselect(df$pov, wmin = .05, wstep = .05, nwindows = 100)

# (iii) continuity-in-density tests
rddensity(df$pov, all = T) %>% summary
rddensity(df$pov, all = T, h = c(5,5)) %>% summary
rddensity(df$pov, all = T, h = c(1,1)) %>% summary

###############################################################################
# Q2.2.1: constant treatment effect model

# generate variables for polynomials of order 3-6
df %<>% mutate(p2=pov^2, p3=pov^3, p4=pov^4, p5=pov^5, p6=pov^6)

# run regressions
r2_1a <- lm(rel_post ~ t + pov + p2 + p3, df)
r2_1b <- lm(rel_post ~ t + pov + p2 + p3 + p4, df)
r2_1c <- lm(rel_post ~ t + pov + p2 + p3 +p4 + p5, df)
r2_1d <- lm(rel_post ~ t + pov + p2 + p3 + p4 + p5 + p6, df)

# grab estimates for table
tab <- matrix(NA,2,4) %>% as.data.frame()
rownames(tab) <- c('Point Estimate', 'Standard Error')
colnames(tab) <- c('p=3', 'p=4', 'p=5', 'p=6')

# point estimates
tab[1,1] <- r2_1a$coefficients['t']
tab[1,2] <- r2_1b$coefficients['t']
tab[1,3] <- r2_1c$coefficients['t']
tab[1,4] <- r2_1d$coefficients['t']

# robust standard errors
tab[2,1] <- diag(vcovHC(r2_1a, type = "HC2")) %>% sqrt() %>% .['t']
tab[2,2] <- diag(vcovHC(r2_1b, type = "HC2")) %>% sqrt() %>% .['t']
tab[2,3] <- diag(vcovHC(r2_1c, type = "HC2")) %>% sqrt() %>% .['t']
tab[2,4] <- diag(vcovHC(r2_1d, type = "HC2")) %>% sqrt() %>% .['t']

# write to LaTeX
xtable(tab)

# generate dataframes for plotting
d2_1a <- data.frame(pov = df$pov, pred = r2_1a$fitted.values)
d2_1b <- data.frame(pov = df$pov, pred = r2_1b$fitted.values)
d2_1c <- data.frame(pov = df$pov, pred = r2_1c$fitted.values)
d2_1d <- data.frame(pov = df$pov, pred = r2_1d$fitted.values)

# generate plot for each order
p2_1a <- ggplot(d2_1a, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 3')
p2_1b <- ggplot(d2_1b, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 4')
p2_1c <- ggplot(d2_1c, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 5')
p2_1d <- ggplot(d2_1d, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 6')

# combine plots
png('r/2_1.png')
multiplot(p2_1a, p2_1c, p2_1b, p2_1d, cols=2)
dev.off()

###############################################################################
# Q2.2.2: heterogeneous treatment effect model

###############################################################################
# Q2.2.3: local parametric model

###############################################################################
# Q2.3.1: MSE-optimal RD estimators

###############################################################################
# Q2.3.2: Robustness checks

###############################################################################
# Q2.4: Local Randomization Methods

###############################################################################