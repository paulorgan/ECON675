###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS1
# last Update: Sept 12, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

require(tidyverse)
require(magrittr)

# for table output
require(xtable)

# permutation tests
require(perm)

# power calculatoins
require(pwr)

select = dplyr::select

setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS1')

###############################################################################
# read in data
df <- read_csv('LaLonde_1986.csv')

###############################################################################
# Question 2: Implementing Least-Squares Estimators

# add necessary variables
df %<>% mutate(educ2 = educ*educ, blackXearn74 = black * earn74)

reg <- lm(earn78 ~ treat + black + age + educ + educ2 +
            earn74 + blackXearn74 + u74 + u75, data = df)

# start table
q2 <- xtable(reg)
colnames(q2) <- c('coefficient', 'std_error', 't_stat', 'p_val')

# add confidence interval
q2$conf_int <-
  paste0('[',round(q2$coefficient-1.96*q2$std_error,2),', ',
         round(q2$coefficient+1.96*q2$std_error,2),']')

# round cols
xtable(q2)

###############################################################################
# Question 3: Analysis of Experiments
# 1) Neyman's Approach:

# 1a) ATE
N1 = sum(df$treat==1)
N0 = sum(df$treat==0)

sumY1 = sum(df$earn78[df$treat==1])
sumY0 = sum(df$earn78[df$treat==0])

Ybar1 = sumY1/N1
Ybar0 = sumY0/N0

ATE = Ybar1-Ybar0

# 1b) t-test
S1 = (1/(N1-1))*var(df$earn78[df$treat==1])
S0 = (1/(N0-1))*var(df$earn78[df$treat==0])

Tstat = ATE / sqrt(S1+S0)
pval = 2*pnorm(-abs(Tstat))

CI = paste0('[',round(ATE-1.96*sqrt(S1+S0),2),', ',round(ATE+1.96*sqrt(S1+S0),2),']')

# Canned version for comparison
ttest <- t.test(earn78 ~ treat, data = df)
ttest

# 2) Fisher's Approach
# Fisher
fisher_1 <- permTS(earn78 ~ treat, data = df,
                   alternative = 'two.sided', method = 'exact.mc',
                   control = permControl(nmc=999,p.conf.level=.95))
fisher_1

# Kolgomorov-Smirnov
earn78_0 <- df$earn78[df$treat==0]
earn78_1 <- df$earn78[df$treat==1]
ks.test(earn78_0, earn78_1, alternative = 'two.sided', exact = T)

# 3) Power Calculations
pwr.t2n.test(N0,N1,d=.5,sig.level=.05)

###############################################################################