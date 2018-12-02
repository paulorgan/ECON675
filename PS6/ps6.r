###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS6
# Last Update: Dec 2, 2018
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

rm(p1_2i); gc()

###############################################################################
# Q2.2.1: constant treatment effect model (additive separable effect)

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

# clean up
rm(list=ls(pattern="r2_|d2_|p2_|tab")); gc()


###############################################################################
# Q2.2.2: heterogeneous treatment effect model (fully interacted effect)

# define new interacted terms for polynomials
df %<>% mutate(p1t = pov*t, p1u = pov*(1-t), p2t = p2*t, p2u = p2*(1-t),
               p3t = p3*t, p3u = p3*(1-t), p4t = p4*t, p4u = p4*(1-t),
               p5t = p5*t, p5u = p5*(1-t), p6t = p6*t, p6u = p6*(1-t))

# define reg equations
v3 <- c('t', 'p1t', 'p1u', 'p2t', 'p2u', 'p3t', 'p3u')
v4 <- c(v3, 'p4t', 'p4u')
v5 <- c(v4, 'p5t', 'p5u')
v6 <- c(v5, 'p6t', 'p6u')

f3 <- paste('rel_post', paste(v3, collapse = ' + '), sep = ' ~ ') %>% as.formula()
f4 <- paste('rel_post', paste(v4, collapse = ' + '), sep = ' ~ ') %>% as.formula()
f5 <- paste('rel_post', paste(v5, collapse = ' + '), sep = ' ~ ') %>% as.formula()
f6 <- paste('rel_post', paste(v6, collapse = ' + '), sep = ' ~ ') %>% as.formula()

# run regressions
r2_2a <- lm(f3, df); r2_2b <- lm(f4, df); r2_2c <- lm(f5, df); r2_2d <- lm(f6, df)

# grab estimates for table
tab <- matrix(NA,2,4) %>% as.data.frame()
rownames(tab) <- c('Point Estimate', 'Standard Error')
colnames(tab) <- c('p=3', 'p=4', 'p=5', 'p=6')

# point estimates
tab[1,1] <- r2_2a$coefficients['t']
tab[1,2] <- r2_2b$coefficients['t']
tab[1,3] <- r2_2c$coefficients['t']
tab[1,4] <- r2_2d$coefficients['t']

# robust standard errors
tab[2,1] <- diag(vcovHC(r2_2a, type = "HC2")) %>% sqrt() %>% .['t']
tab[2,2] <- diag(vcovHC(r2_2b, type = "HC2")) %>% sqrt() %>% .['t']
tab[2,3] <- diag(vcovHC(r2_2c, type = "HC2")) %>% sqrt() %>% .['t']
tab[2,4] <- diag(vcovHC(r2_2d, type = "HC2")) %>% sqrt() %>% .['t']

# write to LaTeX
xtable(tab)

# generate dataframes for plotting
d2_2a <- data.frame(pov = df$pov, pred = r2_2a$fitted.values)
d2_2b <- data.frame(pov = df$pov, pred = r2_2b$fitted.values)
d2_2c <- data.frame(pov = df$pov, pred = r2_2c$fitted.values)
d2_2d <- data.frame(pov = df$pov, pred = r2_2d$fitted.values)

# generate plot for each order
p2_2a <- ggplot(d2_2a, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 3')
p2_2b <- ggplot(d2_2b, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 4')
p2_2c <- ggplot(d2_2c, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 5')
p2_2d <- ggplot(d2_2d, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 6')

# combine plots
png('r/2_2.png')
multiplot(p2_2a, p2_2c, p2_2b, p2_2d, cols=2)
dev.off()

# clean up
rm(list=ls(pattern="r2_|d2_|p2_|v[0-9]|f[0-9]|tab")); gc()

###############################################################################
# Q2.2.3: local parametric model (p = 0,1,2 and h = 1,5,9,18)

# for h = 1
r2_3h1a <- lm(rel_post ~ t, df %>% filter(abs(pov) <= 1))
r2_3h1b <- lm(rel_post ~ t + p1t + p1u, df %>% filter(abs(pov) <= 1))
r2_3h1c <- lm(rel_post ~ t + p1t + p1u + p2t + p2u, df %>% filter(abs(pov) <= 1))

# grab estimates for table
tab <- matrix(NA,12,4) %>% as.data.frame()
tab[,1] <- c('h=1','b','se','h=5','b','se',
             'h=9','b', 'se','h=18', 'b', 'se')
colnames(tab) <- c('','p=0', 'p=1', 'p=2')

# point estimates
tab[2,2] <- r2_3h1a$coefficients['t']
tab[2,3] <- r2_3h1b$coefficients['t']
tab[2,4] <- r2_3h1c$coefficients['t']

# robust standard errors
tab[3,2] <- diag(vcovHC(r2_3h1a, type = "HC2")) %>% sqrt() %>% .['t']
tab[3,3] <- diag(vcovHC(r2_3h1b, type = "HC2")) %>% sqrt() %>% .['t']
tab[3,4] <- diag(vcovHC(r2_3h1c, type = "HC2")) %>% sqrt() %>% .['t']

# generate dataframes for plotting
d2_3h1a <- data.frame(pov = df$pov[abs(df$pov) <= 1], pred = r2_3h1a$fitted.values)
d2_3h1b <- data.frame(pov = df$pov[abs(df$pov) <= 1], pred = r2_3h1b$fitted.values)
d2_3h1c <- data.frame(pov = df$pov[abs(df$pov) <= 1], pred = r2_3h1c$fitted.values)

# generate plot for each order
p2_3h1a <- ggplot(d2_3h1a, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 0, h = 1')
p2_3h1b <- ggplot(d2_3h1b, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 1, h = 1')
p2_3h1c <- ggplot(d2_3h1c, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 2, h = 1')

# combine plots
comb <- grid.arrange(p2_3h1a, p2_3h1b, p2_3h1c, ncol=3)
ggsave('r/2_3h1.png', plot = comb)

# for h = 5
r2_3h5a <- lm(rel_post ~ t, df %>% filter(abs(pov) <= 5))
r2_3h5b <- lm(rel_post ~ t + p1t + p1u, df %>% filter(abs(pov) <= 5))
r2_3h5c <- lm(rel_post ~ t + p1t + p1u + p2t + p2u, df %>% filter(abs(pov) <= 5))

# point estimates
tab[5,2] <- r2_3h5a$coefficients['t']
tab[5,3] <- r2_3h5b$coefficients['t']
tab[5,4] <- r2_3h5c$coefficients['t']

# robust standard errors
tab[6,2] <- diag(vcovHC(r2_3h5a, type = "HC2")) %>% sqrt() %>% .['t']
tab[6,3] <- diag(vcovHC(r2_3h5b, type = "HC2")) %>% sqrt() %>% .['t']
tab[6,4] <- diag(vcovHC(r2_3h5c, type = "HC2")) %>% sqrt() %>% .['t']

# generate dataframes for plotting
d2_3h5a <- data.frame(pov = df$pov[abs(df$pov) <= 5], pred = r2_3h5a$fitted.values)
d2_3h5b <- data.frame(pov = df$pov[abs(df$pov) <= 5], pred = r2_3h5b$fitted.values)
d2_3h5c <- data.frame(pov = df$pov[abs(df$pov) <= 5], pred = r2_3h5c$fitted.values)

# generate plot for each order
p2_3h5a <- ggplot(d2_3h5a, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 0, h = 5')
p2_3h5b <- ggplot(d2_3h5b, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 1, h = 5')
p2_3h5c <- ggplot(d2_3h5c, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 2, h = 5')

# combine plots
comb <- grid.arrange(p2_3h5a, p2_3h5b, p2_3h5c, ncol=3)
ggsave('r/2_3h5.png', plot = comb)

# for h = 9
r2_3h9a <- lm(rel_post ~ t, df %>% filter(abs(pov) <= 9))
r2_3h9b <- lm(rel_post ~ t + p1t + p1u, df %>% filter(abs(pov) <= 9))
r2_3h9c <- lm(rel_post ~ t + p1t + p1u + p2t + p2u, df %>% filter(abs(pov) <= 9))

# point estimates
tab[8,2] <- r2_3h5a$coefficients['t']
tab[8,3] <- r2_3h5b$coefficients['t']
tab[8,4] <- r2_3h5c$coefficients['t']

# robust standard errors
tab[9,2] <- diag(vcovHC(r2_3h9a, type = "HC2")) %>% sqrt() %>% .['t']
tab[9,3] <- diag(vcovHC(r2_3h9b, type = "HC2")) %>% sqrt() %>% .['t']
tab[9,4] <- diag(vcovHC(r2_3h9c, type = "HC2")) %>% sqrt() %>% .['t']

# generate dataframes for plotting
d2_3h9a <- data.frame(pov = df$pov[abs(df$pov) <= 9], pred = r2_3h9a$fitted.values)
d2_3h9b <- data.frame(pov = df$pov[abs(df$pov) <= 9], pred = r2_3h9b$fitted.values)
d2_3h9c <- data.frame(pov = df$pov[abs(df$pov) <= 9], pred = r2_3h9c$fitted.values)

# generate plot for each order
p2_3h9a <- ggplot(d2_3h9a, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 0, h = 9')
p2_3h9b <- ggplot(d2_3h9b, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 1, h = 9')
p2_3h9c <- ggplot(d2_3h9c, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 2, h = 9')

# combine plots
comb <- grid.arrange(p2_3h9a, p2_3h9b, p2_3h9c, ncol=3)
ggsave('r/2_3h9.png', plot = comb)

# for h = 18
r2_3h18a <- lm(rel_post ~ t, df %>% filter(abs(pov) <= 18))
r2_3h18b <- lm(rel_post ~ t + p1t + p1u, df %>% filter(abs(pov) <= 18))
r2_3h18c <- lm(rel_post ~ t + p1t + p1u + p2t + p2u, df %>% filter(abs(pov) <= 18))

# point estimates
tab[11,2] <- r2_3h18a$coefficients['t']
tab[11,3] <- r2_3h18b$coefficients['t']
tab[11,4] <- r2_3h18c$coefficients['t']

# robust standard errors
tab[12,2] <- diag(vcovHC(r2_3h18a, type = "HC2")) %>% sqrt() %>% .['t']
tab[12,3] <- diag(vcovHC(r2_3h18b, type = "HC2")) %>% sqrt() %>% .['t']
tab[12,4] <- diag(vcovHC(r2_3h18c, type = "HC2")) %>% sqrt() %>% .['t']

# generate dataframes for plotting
d2_3h18a <- data.frame(pov = df$pov[abs(df$pov) <= 18], pred = r2_3h18a$fitted.values)
d2_3h18b <- data.frame(pov = df$pov[abs(df$pov) <= 18], pred = r2_3h18b$fitted.values)
d2_3h18c <- data.frame(pov = df$pov[abs(df$pov) <= 18], pred = r2_3h18c$fitted.values)

# generate plot for each order
p2_3h18a <- ggplot(d2_3h18a, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 0, h = 18')
p2_3h18b <- ggplot(d2_3h18b, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 1, h = 18')
p2_3h18c <- ggplot(d2_3h18c, aes(x = pov, y = pred)) + geom_point() +
  theme_minimal() + labs(x = 'Pov Rate 1960 (Normalized)', y = 'Fitted Values',
                         title = 'Order 2, h = 18')

# combine plots
comb <- grid.arrange(p2_3h18a, p2_3h18b, p2_3h18c, ncol=3)
ggsave('r/2_3h18.png', plot = comb)

# combine all and check
comb <- grid.arrange(p2_3h1a, p2_3h1b, p2_3h1c,
                     p2_3h5a, p2_3h5b, p2_3h5c,
                     p2_3h9a, p2_3h9b, p2_3h9c,
                     p2_3h18a, p2_3h18b, p2_3h18c, ncol=3)
ggsave('r/2_3.png', plot = comb)

# write table to LaTeX
xtable(tab)

# clean up
rm(list=ls(pattern="d2_|p2_|r2_|comb|tab")); gc()

###############################################################################
# Q2.3.1: MSE-optimal RD estimators

# just manually typings these into LaTeX
# 'all=T' gives the three methods
# p and q tell the polynomial orders to use
rdrobust(df$rel_post, df$pov, p = 0, q = 1, c = 0, all = T) %>% summary()
rdrobust(df$rel_post, df$pov, p = 1, q = 2, c = 0, all = T) %>% summary()
rdrobust(df$rel_post, df$pov, p = 2, q = 3, c = 0, all = T) %>% summary()

###############################################################################
# Q2.3.2: Robustness checks

# a) Placebo outcome tests

# check preintervention related mortality
rdrobust(df$rel_pre, df$pov, p =1, q = 2, c = 0, all = T) %>% summary()
# check postintervention unrelated mortality
rdrobust(df$inj_post, df$pov, p =1, q = 2, c = 0, all = T) %>% summary()

# b) Bandwidth and Kernel sensitivity
# empty matrix to fill with results
tbl <- matrix(NA,3,10)

# loop over 10 bandwidths, for three kernels, grab point estimates
for(h in 1:10){
  tbl[1,h] <- rdrobust(df$rel_pre, df$pov, p =1, q = 2, c = 0,
                       all = T, kernel = 'tri')$Estimate[2]
  tbl[2,h] <- rdrobust(df$rel_pre, df$pov, p =1, q = 2, c = 0,
                       all = T, kernel = 'uni')$Estimate[2]
  tbl[3,h] <- rdrobust(df$rel_pre, df$pov, p =1, q = 2, c = 0,
                       all = T, kernel = 'epa')$Estimate[2]
}

# clean up and write to LaTeX
tbl %<>% as.data.frame()
rownames(tbl) <- c('Triangular', 'Uniform', 'Epanechnikov')
colnames(tbl) <- 1:10
xtable(tbl)

# c) "donut hole" approach
# sort by proximity to 0, and number them
df %<>% mutate(pov_abs = abs(pov)) %>% arrange(pov_abs)
df$abs_pov_rank <- 1:nrow(df)

# empty matrix to fill with results
tbl <- matrix(NA, 1, 10)

# loop over 10 l's
for(l in 1:10){
  tbl[1,l] <- rdrobust(df$rel_pre[df$abs_pov_rank > l], df$pov[df$abs_pov_rank > l],
                       p = 1, q = 2, c = 0, all = T)$Estimate[2]
}

# clean up and write to LaTeX
tbl %<>% as.data.frame()
colnames(tbl) <- 1:10
xtable(tbl)

# d) Placebo cutoff approach
cutoffs <- c(-10, -8, -6, -4, -2, 2, 4, 6, 8, 10)

# empty matrix to fill with results
tbl <- matrix(NA, 2, 10)

# loop over 10 cutoffs, save point estimates and p-values
for(i in 1:10){
  temp <- rdrobust(df$rel_pre, df$pov, p = 1, q = 2, c = cutoffs[i], all = T)
  tbl[1,i] <- temp$Estimate[2]
  tbl[2,i] <- temp$pv[2]
}

# clean up and write to LaTeX
tbl %<>% as.data.frame()
rownames(tbl) <- c('Point estimate', 'p-value')
colnames(tbl) <- 1:10
xtable(tbl)

# clean up Q2.3
rm(r3_1a, r3_1b, r3_1c, tbl, temp, cutoffs, h, i, l); gc()

###############################################################################
# Q2.4: Local Randomization Methods

###############################################################################