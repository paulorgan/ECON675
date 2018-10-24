###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS4
# Last Update: Oct 24, 2018
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
df %<>% mutate(treat_full = ifelse(treat == 1, 1, 0))

# add other variable we will use in the analysis
df %<>% mutate(logre74 = log(re74 + 1), logre75 = log(re75 + 1),
               age2 = age^2, educ2 = educ^2, age3 = age^3,
               bXu74 = black * u74, eXre74 = educ * logre74)

# empty table to fill with results (1 + 6*3)
ate <- matrix(NA, nrow=19,ncol=8) %>% as.data.frame

# add headings
names(ate) <- c('e_tau', 'e_se', 'e_CI_l', 'e_CI_h', 'p_tau', 'p_se', 'p_CI_l', 'p_CI_h')

# add empty table for att
att <- ate

# define three model specifications
va <- c('age', 'educ', 'black', 'hisp', 'married',
        'nodegr', 'logre74', 'logre75')
vb <- c(va, 'age2', 'educ2', 'u74', 'u75')
vc <- c(vb, 'age3', 'bXu74', 'eXre74')

a <- paste('re78', paste(c('treat_full', va),
                         collapse = ' + '), sep = ' ~ ') %>% as.formula()
b <- paste('re78', paste(c('treat_full', vb),
                         collapse = ' + '), sep = ' ~ ') %>% as.formula()
c <- paste('re78', paste(c('treat_full', vc),
                         collapse = ' + '), sep = ' ~ ') %>% as.formula()

ta <- paste('treat_full', paste(va,collapse = ' + '), sep = ' ~ ') %>% as.formula()
ta <- paste('treat_full', paste(vb, collapse = ' + '), sep = ' ~ ') %>% as.formula()
ta <- paste('treat_full', paste(vc, collapse = ' + '), sep = ' ~ ') %>% as.formula()
  
rm(va, vb, vc)

###############################################################################
# Q2.1: Difference-in-Means

# difference in means for experimental data
r1e <- lm(re78 ~ treat_full, data = df %>% filter(treat != 2))

# extract tau
ate$e_tau[1]  <- r1e$coefficients['treat_full']

# robust standard errors (HC1 is equivalent to Stata robust)
ate$e_se[1]   <- diag(vcovHC(r1e, type = "HC1")) %>% sqrt() %>% .['treat_full']

# difference in means for PSID data
r1p <- lm(re78 ~ treat_full, data = df %>% filter(treat != 0))

# extract tau
ate$p_tau[1]  <- r1p$coefficients['treat_full']

# robust standard errors
ate$p_se[1]   <- diag(vcovHC(r1p, type = "HC1")) %>% sqrt() %>% .['treat_full']

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
ate$e_tau[2]  <- r2ea$coefficients['treat_full']
ate$e_tau[3]  <- r2eb$coefficients['treat_full']
ate$e_tau[4]  <- r2ec$coefficients['treat_full']

# robust standard errors 
ate$e_se[2]   <- diag(vcovHC(r2ea, type = "HC1")) %>% sqrt() %>% .['treat_full']
ate$e_se[3]   <- diag(vcovHC(r2eb, type = "HC1")) %>% sqrt() %>% .['treat_full']
ate$e_se[4]   <- diag(vcovHC(r2ec, type = "HC1")) %>% sqrt() %>% .['treat_full']

# OLS for PSID data
r2pa <- lm(a, data = df %>% filter(treat != 0))
r2pb <- lm(b, data = df %>% filter(treat != 0))
r2pc <- lm(c, data = df %>% filter(treat != 0))

# extract tau
ate$p_tau[2]  <- r2pa$coefficients['treat_full']
ate$p_tau[3]  <- r2pb$coefficients['treat_full']
ate$p_tau[4]  <- r2pc$coefficients['treat_full']

# robust standard errors 
ate$p_se[2]   <- diag(vcovHC(r2pa, type = "HC1")) %>% sqrt() %>% .['treat_full']
ate$p_se[3]   <- diag(vcovHC(r2pb, type = "HC1")) %>% sqrt() %>% .['treat_full']
ate$p_se[4]   <- diag(vcovHC(r2pc, type = "HC1")) %>% sqrt() %>% .['treat_full']

# again, ATE is same as ATT here
att[2:4,] <- ate[2:4,]

###############################################################################
# Q2.3: Regression Imputation


###############################################################################
# Q2.4: Inverse Probability Weighting

###############################################################################
# Q2.5: Doubly Robust

###############################################################################
# Q2.6: Nearest Neighbor Matching

# this may be useful:
# https://cran.r-project.org/web/packages/MatchIt/vignettes/matchit.pdf

# also this:
# https://www.r-bloggers.com/using-the-r-matchit-package-for-propensity-score-analysis/

# don't think this is right...
r6ea_p <- matchit(ta, data = df %>% filter(treat != 2), method = 'nearest')
df_r6ea <- match.data(r6ea_p)
mean(df_r6ea$re78[df_r6ea$treat_full==1])-mean(df_r6ea$re78[df_r6ea$treat_full==0])
t.test(df_r6ea$re78[df_r6ea$treat_full==1],df_r6ea$re78[df_r6ea$treat_full==0])

###############################################################################
# Q2.7: Propensity Score Matching

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

###############################################################################
### Question 3: Post-model Selection Inference
###############################################################################
# Q3 setup

###############################################################################
# Q3.1: Summary Statistics and Kernel Density Plots

###############################################################################
# Q3.2: Empirical Coverage Rates

###############################################################################