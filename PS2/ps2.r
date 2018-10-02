###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS2
# Last Update: Oct 2, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse) # data cleaning and manipulation
require(magrittr)  # syntax
require(ggplot2)   # plots
require(kedd)      # kernel estimation
require(ks)        # kernel estimation
require(car)       # heteroskedastic robust SEs

setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS2')

###############################################################################
## Question 1: Kernel Density Estimation
# Q1.3a

# sample size
n     <- 1000

# equally weight two distributions
set.seed(22)
comps <- sample(1:2,prob=c(.5,.5),size=n,replace=T)

# Normal density specs
mus <- c(-1.5, 1)
sds <- sqrt(c(1.5, 1))

# generate sample
set.seed(22)
samp <- rnorm(n=n,mean=mus[comps],sd=sds[comps])

# check plot
plot(density(samp))

# define Kernel function: K(u)=.75(1-u^2)(ind(abs(u)<=1))
K0 <- function(u){
  out <- .75 * (1-u^2) * (abs(u) <= 1)
}

# first derivative of Kernel function
K1 <- function(u){
  out <- .75 * (-2 * u) * (abs(u) <= 1)
}

# compute optimal bandwidth using kedd package
h.amise(samp, kernel = 'epanechnikov', deriv.order=0)
plot.h.amise(h.amise(samp, kernel = 'epanechnikov', deriv.order= 0))

h.amise(samp, kernel = 'epanechnikov', deriv.order=1)

# what about the np package? see function npregbw
# https://cran.r-project.org/web/packages/np/np.pdf

# trying using ks package
# can't implement epanechnikov kernel here...only Normal
hamise.mixt(mus,sds,props=c(.5,.5),samp=n,deriv.order=0)
hamise.mixt(mus,sds,props=c(.5,.5),samp=n,deriv.order=1)

# Q1.3b
# simulate 1000 times
M <- 1000

###############################################################################
## Question 2: Linear Smoothers, Cross-Validation, and Series
# clean up
rm(list = ls())
gc()

###############################################################################
## Q5a
# define datagenerating process
set.seed(22)
dgp <- function(n){
  # input:  sample size n
  # output: n draws of X and Y according to DGPs as specified in PSet
  
  # X ~ Uniform(-1,1)
  X <- runif(n,-1,1)
  
  # Epsilon ~ x^2 * (Chisq_5 - 5)
  E <- (X^2) * (rchisq(n,5)-5)
  
  # Y ~ exp(-.1*(4x-1)^2) * sin(5x) + Eps
  Y <- exp(-0.1 * (4*X-1)^2) * sin(5*X) + E
  
  out = data.frame(X=X,Y=Y)
}

# sample size
n <- 1000

# replications
M <- 1000

###############################################################################
## Q5b
# series truncation
K  <- 1:20
nK <- length(K)

# define series cross-validation function
series_cv <- function(n, X, Y, nK, K){
  # input:  n draws of X and Y, series truncation K of length nK
  # output: nK prediction errors
  
  # start with polynomial basis of X
  X_poly <- cbind(rep(1,n), poly(X,degree = nK))
  
  # QR decomposition
  X_qr <- qr(X_poly)
  
  # coefficients
  coefs <- qr.coef(X_qr,Y)
  
  # cycle up through each power and store prediction errors
  out <- rep(NA,nK)
  for(k in 1:nK){
    X_poly_k <- X_poly[,1:(K[k]+1)]
    coefs_k  <- matrix(coefs[1:(K[k]+1)],nrow=K[k]+1)
    Y_hat_k  <- X_poly_k %*% coefs_k
    
    # w (see part 1 of question 2)
    w_k <- diag(X_poly_k %*% solve(t(X_poly_k) %*% X_poly_k) %*% t(X_poly_k))
    
    # prediction error
    cv_k <- mean( ((Y-Y_hat_k)/(1-w_k))^2 )
    
    # save
    out[k] <- cv_k
  }
  
  return(out)
}

# run series cv formula 1000 times, save results for plotting
# we want to capture the MSE for each K for each rep
MSEs <- matrix(NA, ncol=nK, nrow=M)
set.seed(22)
for(m in 1:M){
  df <- dgp(n)
  MSEs[m,] <- series_cv(n=n, X=df$X, Y=df$Y, nK=nK, K=K)
}

# average CV(K) across simulations:
averages <- data.frame(K = 1:20, avg_MSE = colMeans(MSEs))

# identify optimal CV estimator
averages %<>% mutate(Group = (avg_MSE == min(avg_MSE)))
K_CV <- averages %>% filter(Group) %>% select(K) %>% as.numeric

# plot average CV(K), highlight the optimal one
plot <- ggplot(data = averages, aes(x = K, y = avg_MSE, color = Group)) +
  geom_point(size = 2) + scale_color_manual(values=c('black','red')) +
  theme(legend.position='none')
plot
ggsave('q5b_R.png')

###############################################################################
## Q5c
# define grid of evaluation points
grid     <- seq(-1,1,.1)
eval_pts <- length(grid)

# we know the true value of the regression function (from dgp defined above)
f_true <- exp(-0.1 * (4*grid-1)^2) * sin(5*grid)

# estimate regression function using polynomial basis (7th degree based on (b))
# generate polynomial basis for grid points
grid_poly <- cbind(1, grid, grid^2, grid^3, grid^4, grid^5, grid^6, grid^7)

# define empty matrices to fill with estimates
ests <- matrix(NA, ncol=eval_pts, nrow=M)
SEs  <- matrix(NA, ncol=eval_pts, nrow=M)

# simulate M times
set.seed(22)
for(m in 1:M){
  # draw data
  df <- dgp(n)
  X <- df$X
  Y <- df$Y
  
  # create polynomial
  X_poly <- cbind(1,X,X^2,X^3,X^4,X^5,X^6,X^7)
  
  # run regression using drawn data
  reg   <- lm(Y ~ X_poly - 1)
  betas <- coefficients(reg)
  vars  <- hccm(reg, type='hc0') # heteroskedasticity-corrected using 'car' package
  
  # calculate and save estimates and SEs using grid evaluation points
  ests[m,] <- grid_poly %*% betas
  SEs[m,]  <- sqrt(diag(grid_poly %*% vars %*% t(grid_poly)))
}

# average across the simulations to create estimated regression function
f_est <- colMeans(ests)

# calculate average confidence intervals across the simulations
f_est_low  <- colMeans(ests)-1.96*colMeans(SEs)
f_est_high <- colMeans(ests)+1.96*colMeans(SEs)

# gather data for plotting
df <- data.frame(x = grid, true = f_true, est = f_est,
                 CI_low = f_est_low, CI_high = f_est_high) %>%
  gather(key = type, value = value, -x)

# plot in one graph, note specification of options is alphabetical by type
# types = CI_high, CI_low, est, true
plot <- ggplot(data = df, aes(x = x, y = value, color = type)) +
  geom_line(aes(linetype=type, color=type)) +
  geom_point(aes(color=type, size=type)) +
  scale_linetype_manual(values = c('dotted', 'dotted', 'longdash', 'solid')) +
  scale_color_manual(values = c('black', 'black', 'red', 'blue')) +
  scale_size_manual(values = c(0,0,2,2))
plot
ggsave('q5c_R.png')

###############################################################################
## Q5d
# Now estimating the derivative of the regression function
# This will reuse my code from above, with new polynomials (taking derivs)

#  true value of derivative of regression function (product rule)
f1_true <- exp(-.1*(4*grid-1)^2)*5*cos(5*grid) + sin(5*grid)*(-.8*(4*grid-1)*exp(-.1*(4*grid-1)^2))

# estimate regression function using polynomial basis (7th degree based on (b))
# generate polynomial basis for grid points, first derivatives of grid_poly
grid1_poly <- cbind(0, 1, 2*grid, 3*grid^2, 4*grid^3, 5*grid^4, 6*grid^5, 7*grid^6)

# define empty matrices to fill with estimates
ests1 <- matrix(NA, ncol=eval_pts, nrow=M)
SEs1  <- matrix(NA, ncol=eval_pts, nrow=M)

# simulate M times
set.seed(22)
for(m in 1:M){
  # draw data
  df <- dgp(n)
  X <- df$X
  Y <- df$Y
  
  # create polynomial
  X_poly <- cbind(1,X,X^2,X^3,X^4,X^5,X^6,X^7)
  
  # run regression using drawn data
  reg   <- lm(Y ~ X_poly - 1)
  betas <- coefficients(reg)
  vars  <- hccm(reg, type='hc0') # heteroskedasticity-corrected using 'car' package
  
  # calculate and save estimates and SEs using grid evaluation points
  ests1[m,] <- grid1_poly %*% betas
  SEs1[m,]  <- sqrt(diag(grid1_poly %*% vars %*% t(grid1_poly)))
}

# average across the simulations to create estimated regression function
f1_est <- colMeans(ests1)

# calculate average confidence intervals across the simulations
f1_est_low  <- colMeans(ests1)-1.96*colMeans(SEs1)
f1_est_high <- colMeans(ests1)+1.96*colMeans(SEs1)

# gather data for plotting
df1 <- data.frame(x = grid, true = f1_true, est = f1_est,
                  CI_low = f1_est_low, CI_high = f1_est_high) %>%
   gather(key = type, value = value, -x)

# plot in one graph, note specification of options is alphabetical by type
# types = CI_high, CI_low, est, true
plot <- ggplot(data = df1, aes(x = x, y = value, color = type)) +
  geom_line(aes(linetype=type, color=type)) +
  geom_point(aes(color=type, size=type)) +
  scale_linetype_manual(values = c('dotted', 'dotted', 'longdash', 'solid')) +
  scale_color_manual(values = c('black', 'black', 'red', 'blue')) +
  scale_size_manual(values = c(0,0,2,2))
plot
ggsave('q5d_R.png')

###############################################################################