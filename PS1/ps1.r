###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS1
# last Update: Sept 17, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse)
require(magrittr)
require(xtable) # for table output
require(perm) # permutation tests
require(boot) # bootstrapping
require(ggplot2) # plots

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
## Question 3: Analysis of Experiments
# 3.1) Neyman's Approach:

# 3.1a) ATE
N1 = sum(df$treat==1)
N0 = sum(df$treat==0)

sumY1 = sum(df$earn78[df$treat==1])
sumY0 = sum(df$earn78[df$treat==0])

Ybar1 = sumY1/N1
Ybar0 = sumY0/N0

ATE = Ybar1-Ybar0

# 3.1b) t-test
S1 = (1/(N1-1))*var(df$earn78[df$treat==1])
S0 = (1/(N0-1))*var(df$earn78[df$treat==0])

se <- sqrt(S1+S0)

Tstat = ATE / se
pval = 2*pnorm(-abs(Tstat))

CI_31b = paste0('[',round(ATE-1.96*sqrt(S1+S0),2),', '
            ,round(ATE+1.96*sqrt(S1+S0),2),']')

# Canned version for comparison
ttest <- t.test(earn78 ~ treat, data = df)
ttest

###############################################################################
# 3.2) Fisher's Approach
# 3.2a) p-Value
# Fisher
fisher_1 <- permTS(earn78 ~ treat, data = df,
                   alternative = 'two.sided', method = 'exact.mc',
                   control = permControl(nmc=999,p.conf.level=.95))
fisher_1

# Kolgomorov-Smirnov
earn78_0 <- df$earn78[df$treat==0]
earn78_1 <- df$earn78[df$treat==1]
ks.test(earn78_0, earn78_1, alternative = 'two.sided', exact = T)

# 3.2b) Confidence Interval
# Imputation assuming ATE is constant
# (Generating Yi(1) and Yi(0) for each i, assuming ATE estimate is constant)
Y1_imp <- (df$treat==1) * df$earn78 + (df$treat==0) * (df$earn78 + ATE)
Y0_imp <- (df$treat==1) * (df$earn78 - ATE) + (df$treat==0) * df$earn78

# define statistic (difference in means)
boot_T <- function(x, ind) {
  mean(Y1_imp[df$treat[ind]==1]) - mean(Y0_imp[df$treat[ind]==0])
}

# run bootstrap using defined statistic
boot_results <- boot(df, R = 999, statistic = boot_T,
                     sim = "permutation", stype = "i")

# construct 95% confidence interval
CI_32b = paste0('[', quantile(boot_results$t,0.025), ', ',
               quantile(boot_results$t,0.975), ']')

###############################################################################
# 3.3) Power Calculations
# 3.3a) graphing power function

# Z value for 95%
Z <- 1.96

# testing tau_0 = 0. plot tau on either side of 0
df_p <- data.frame(tau = seq(-2500, 2500, 25))

# calculate propability of rejection under each alternative tau
df_p %<>% mutate(prob_rej = pnorm(Z-tau/se, lower.tail=F) +
                   pnorm(Z+tau/se, lower.tail=F))

plot <- ggplot(df_p, aes(x = tau, y = prob_rej)) +
  geom_point() + geom_smooth() +
  ylab('Power') + xlab('tau') +
  geom_hline(yintercept=0.05,linetype='dashed',color='red')
ggsave('power_R.png', plot)

# 3.3b) determining minimum sample size
df_n <- data.frame(n = seq(100,5000,5))

# given: probablity of treatment is 2/3
p <- 2/3

# effect we want to detect
tau_0 <- 1000

# variances of observed data
V1 <- var(df$earn78[df$treat==1])
V0 <- var(df$earn78[df$treat==0])

# simulate with different sample sizes
df_n %<>% mutate(n1 = n*p,
                 n0 = n-n1,
                 std_err = sqrt(V1/n1 + V0/n0),
                 power = pnorm(Z-tau_0/std_err, lower.tail=F) +
                   pnorm(Z+tau_0/std_err, lower.tail=F)) 

# find sample size st power is at least .8
min_n <- min(df_n$n[df_n$power>=0.8])

###############################################################################