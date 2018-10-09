###############################################################################
# Author: Paul R. Organ
# Purpose: ECON 675, PS2
# Last Update: Oct 9, 2018
###############################################################################
# Preliminaries
options(stringsAsFactors = F)

# packages
require(tidyverse) # data cleaning and manipulation
require(magrittr)  # syntax
require(ggplot2)   # plots
require(kedd)      # kernel bandwidth estimation
require(car)       # heteroskedastic robust SEs
require(xtable)    # tables for LaTeX

setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS2')

###############################################################################
## Question 1: Kernel Density Estimation
## Q1.3a

# sample size
n     <- 1000

# data generating process
dgp <- function(n){
  # equally weight two distributions
  comps <- sample(1:2,prob=c(.5,.5),size=n,replace=T)
  
  # Normal density specs
  mus <- c(-1.5, 1)
  sds <- sqrt(c(1.5, 1))
  
  # generate sample
  samp <- rnorm(n=n,mean=mus[comps],sd=sds[comps])
  
  return(samp)
}

# true dgp
f_true <- function(x){.5*dnorm(x,-1.5,1.5)+.5*dnorm(x,1,1)}

# second deriv of normal dist
norm_2d <- function(u,meanu,sdu){dnorm(u,mean=meanu,sd=sdu)*
                    (((u-meanu)^2/(sdu^4))-(1/(sdu^2)))}

# f for integration (theoretical)
f_int <- function(x){
  return( (.5*norm_2d(x,-1.5,sqrt(1.5)) + .5*norm_2d(x,1,1))^2 )
}

# function to calculate (theoretically or empirical) optimal h
optimal_h <- function(x,mu,sd){
  k1 <- .75^2 * (2- 4/3 + 2/5)
  k2 <- .75 * (2/3 - 2/5)

  if(x=='theoretical'){
    f <- f_int
    k3 <- integrate(f,-Inf,Inf)$val
  } else{
    f <- function(x,mu,sd){norm_2d(x,mu,sd)^2}
    k3 <- integrate(f,-Inf,Inf,mu=mu,sd=sd)$val
  }
    
  h <- (k1/(k3*k2^2)*(1/n))^(1/5)
  return(h)
}

# theoretically optimal h
h_aimse <- optimal_h('theoretical',NA,NA)

###############################################################################
## Q1.3b
# define Kernel function: K(u)=.75(1-u^2)(ind(abs(u)<=1))
K0 <- function(u){
  out <- .75 * (1-u^2) * (abs(u) <= 1)
}

# function to calculate IMSE
imse <- function(h,X){
  # empty vectors to fill with results
  e_li <- rep(NA,n)
  e_lo <- rep(NA,n)
  
  # loop over each i to do leave one out
  for(i in 1:n){
    # repeat observation for each simulation
    Xi_n <- rep(X[i], n)
    # apply kernel function to (x-x_i)/h
    df <- K0((Xi_n-X)/h)
    
    # fhat with i in
    fhat_li <- mean(df)/h
    # fhat with i out
    fhat_lo <- mean(df[-i])/h
    
    # f(x_i)
    f_xi <- f_true(X[i])
    
    # mse with i in
    e_li[i] <- (fhat_li-f_xi)^2
    # mse with i out
    e_lo[i] <- (fhat_lo-f_xi)^2
  }
  out <- c(mean(e_li),mean(e_lo))
  return(out)
}

# simulate 1000 times
M <- 1000

# sequence of h's to test
hs <- seq(.5,1.5,.1) * h_aimse
nh <- length(hs)

# empty matrices to fill
imse_li <- matrix(NA, nrow=M, ncol=nh)
imse_lo <- matrix(NA, nrow=M, ncol=nh)

# generate matrix with M rows of sampled data
set.seed(22)
for(m in 1:M){
  X <- dgp(n)
  for(j in 1:nh){
    temp <- imse(hs[j],X)
    imse_li[m,j] <- temp[1]
    imse_lo[m,j] <- temp[2]
  }
}

df <- data.frame(h = hs, leavein = colMeans(imse_li), leaveout = colMeans(imse_lo)) %>%
  gather(key = inout, value = imse, -h)

# plot
p <- ggplot(df, aes(x=h,y=imse,color=inout)) + geom_smooth(se=F) +
  theme_minimal()
p
ggsave('q1_3b_R.png')

# averages
h_hat_li <- df %>% filter(inout == 'leavein') %>%
  filter(imse == min(imse)) %>% select(h) %>% as.numeric
h_hat_lo <- df %>% filter(inout == 'leaveout') %>%
  filter(imse == min(imse)) %>% select(h) %>% as.numeric

###############################################################################
## Q1.3d

opt_hs <- rep(NA,M)

set.seed(22)
for(m in 1:M){
  X <- dgp(n)
  mu <- mean(X)
  sd <- sd(X)
  opt_hs[m] <- optimal_h(n,mu,sd)
}

hbar <- mean(opt_hs)

###############################################################################
## Question 2: Linear Smoothers, Cross-Validation, and Series
# clean up
rm(list = ls())
gc()

###############################################################################
## Q2.5a
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
## Q2.5b
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
ggsave('q2_5b_R.png')

###############################################################################
## Q2.5c
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
ggsave('q2_5c_R.png')

###############################################################################
## Q2.5d
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
ggsave('q2_5d_R.png')

###############################################################################
## Question 3: Semiparametric Semi-Linear Model
# clean up
rm(list = ls())
gc()

###############################################################################
## Q3.4a
# define data generating process

# sample size
n <- 500
# replications
M <- 1000

# define data generating process
dgp <- function(n){
  # input:  sample size n
  # output: n draws of X and Y according to DGPs as specified in PSet
  
  # X is a d(5) by n matrix ~ U(-1,1)
  X <- matrix(runif(n*5,-1,1), ncol=5)
  
  # V ~ N(0,1) and U ~ N(0,1)
  V <- rnorm(n)
  U <- rnorm(n)
  
  # Eps = .36...*(1+||X||^2)*V
  E <- 0.3637899*(1+diag(X %*% t(X)))*V
  
  # g_0(X) = exp(||X||^2)
  G <- exp(diag(X %*% t(X)))
  
  # T = ind(||x||+u >= 0) (times 1 to convert from Boolean to numeric)
  Tee <- matrix((sqrt(diag(X %*% t(X))) + U >= 0)*1, ncol = 1)
  
  # Y as defined in problem (assuming \theta_0 = 1)
  Y <- matrix(Tee + G + E, ncol = 1)
  
  # returning list with matrix X, vector Y, vector T
  out = list(X=X,Y=Y,Tee=Tee)
}

# define polynomial basis
K <- c(6,11,21,26,56,61,126,131,252,257,262,267,272,277)

polybasis <- function(X,K){
  # inputs: data matrix X, 'order' K
  # outputs: polynomial basis of 'order' K
  # Note: this gets really ugly at the end,
  # but I was tired and gave up trying to find a more elegant way
  if(K==6){basis <- poly(X,degree=1,raw=T)}
  if(K==11){basis <- cbind(poly(X,degree=1,raw=T),
                           X[,1]^2,X[,2]^2,X[,3]^2,X[,4]^2,X[,5]^2)}
  if(K==21){basis <- poly(X,degree=2,raw=T)}
  if(K==26){basis <- cbind(poly(X,degree=2,raw=T),
                           X[,1]^3,X[,2]^3,X[,3]^3,X[,4]^3,X[,5]^3)}
  if(K==56){basis <- poly(X,degree=3,raw=T)}
  if(K==61){basis <- cbind(poly(X,degree=3,raw=T),
                           X[,1]^4,X[,2]^4,X[,3]^4,X[,4]^4,X[,5]^4)}
  if(K==126){basis <- poly(X,degree=4,raw=T)}
  if(K==131){basis <- cbind(poly(X,degree=4,raw=T),
                            X[,1]^5,X[,2]^5,X[,3]^5,X[,4]^5,X[,5]^5)}
  if(K==252){basis <- poly(X,degree=5,raw=T)}
  if(K==257){basis <- cbind(poly(X,degree=5,raw=T),
                            X[,1]^6,X[,2]^6,X[,3]^6,X[,4]^6,X[,5]^6)}
  if(K==262){basis <- cbind(poly(X,degree=5,raw=T),
                            X[,1]^6,X[,2]^6,X[,3]^6,X[,4]^6,X[,5]^6,
                            X[,1]^7,X[,2]^7,X[,3]^7,X[,4]^7,X[,5]^7)}
  if(K==267){basis <- cbind(poly(X,degree=5,raw=T),
                            X[,1]^6,X[,2]^6,X[,3]^6,X[,4]^6,X[,5]^6,
                            X[,1]^7,X[,2]^7,X[,3]^7,X[,4]^7,X[,5]^7,
                            X[,1]^8,X[,2]^8,X[,3]^8,X[,4]^8,X[,5]^8)}
  if(K==272){basis <- cbind(poly(X,degree=5,raw=T),
                            X[,1]^6,X[,2]^6,X[,3]^6,X[,4]^6,X[,5]^6,
                            X[,1]^7,X[,2]^7,X[,3]^7,X[,4]^7,X[,5]^7,
                            X[,1]^8,X[,2]^8,X[,3]^8,X[,4]^8,X[,5]^8,
                            X[,1]^9,X[,2]^9,X[,3]^9,X[,4]^9,X[,5]^9)}
  if(K==277){basis <- cbind(poly(X,degree=5,raw=T),
                            X[,1]^6,X[,2]^6,X[,3]^6,X[,4]^6,X[,5]^6,
                            X[,1]^7,X[,2]^7,X[,3]^7,X[,4]^7,X[,5]^7,
                            X[,1]^8,X[,2]^8,X[,3]^8,X[,4]^8,X[,5]^8,
                            X[,1]^9,X[,2]^9,X[,3]^9,X[,4]^9,X[,5]^9,
                            X[,1]^10,X[,2]^10,X[,3]^10,X[,4]^10,X[,5]^10)}
  return(basis)
}

###############################################################################
## Q3.4b

# number of different orders to test
nK <- length(K)

# define blank matrices to fill with simulated results
thetas <- matrix(NA, ncol = nK, nrow = M)
SEs    <- matrix(NA, ncol = nK, nrow = M)

set.seed(22)
ptm <- proc.time()
for(m in 1:M){
  # draw data
  data <- dgp(n)
  X   <- data$X
  Y   <- data$Y
  Tee <- data$Tee
  
  # cycle through K orders
  for(k in 1:nK){
    # generate basis, add intercept
    X_poly <- cbind(1,polybasis(X,K[k]))

    # define M_P (I-P(P'P)^{-1}P)
    M_P <- diag(n) - (X_poly %*% solve((t(X_poly) %*% X_poly)) %*% t(X_poly))
    
    # estimate theta(K)
    theta <- (t(Tee) %*% M_P %*% Y)/(t(Tee) %*% M_P %*% Tee)
    
    # sigma (for variance estimate)
    sigma <- diag( as.numeric((M_P %*% (Y - Tee*as.numeric(theta))))^2 )
    
    # standard error
    bread <- solve((t(Tee) %*% M_P %*% Tee))
    se <- sqrt( bread %*% (t(Tee)%*%M_P%*%sigma%*%M_P%*%Tee) %*% bread)
    
    # save to matrix
    thetas[m,k] <- theta
    SEs[m,k]    <- se
  }
}
proc.time() - ptm
# 12 minute runtime

# calculate averages, variances, etc. of simulated values
summ <- matrix(NA, ncol=6, nrow = nK)
for(k in 1:nK){
  summ[k,1] <- K[k] # 'order' K
  summ[k,2] <- mean(thetas[,k]) # avergage theta(K)
  summ[k,3] <- summ[k,2]-1 # average bias of theta(K) (assuming theta_0=1)
  summ[k,4] <- sd(thetas[,k])^2 # sample variance of theta(K)
  summ[k,5] <- mean( (SEs[,k])^2 ) # average of vhat
  # check if CIs for each simulation include 1 (here checking the boundaries)
  summ[k,6] <- 1 - mean( thetas[,k]-1.96*SEs[k] > 1 | thetas[,k]+1.96*SEs[k] < 1)
}

# format
summ %<>% as.data.frame
names(summ) <- c('K', 'avg_theta', 'avg_bias',
                 'samp_variance', 'avg_vhat', 'coverage_rate')

# write for inclusion in latex document
print(xtable(summ,digits=c(0,0,3,3,3,3,3)),include.rownames=F)

###############################################################################
## Q3.4c

# define crossvalidation function
crossval <- function(X, Y, Tee, nK, K){
  # blank vector to fill with MSE
  MSEs <- rep(NA, nK)
  
  # loop through each K to identify optimal bandwidth
  for(k in 1:nK){
    # define polynomial basis
    X_poly <- cbind(1, Tee, polybasis(X, K[k]))
    # QR decomposition
    X_poly_Q <- qr.Q(qr(X_poly))
    XX <- X_poly_Q %*% t(X_poly_Q)
    Y_hat <- XX %*% Y
    W <- diag(XX)
    MSEs[k] <- mean( ((Y-Y_hat)/(1-W))^2 )
  }

  # return the optimal K
  return(K[which.min(MSEs)])
}

# define blank vectors to fill with simulated results and optimal Ks
# (not matrices now, since we are using optimal K)
thetas <- rep(NA, M)
SEs    <- rep(NA, M)
Ks     <- rep(NA, M)

# simulate M times
set.seed(22)
ptm <- proc.time()
for(m in 1:M){
  # draw data
  data <- dgp(n)
  X   <- data$X
  Y   <- data$Y
  Tee <- data$Tee
  
  # given data, estimate optimal K using cross validation
  K_CV <- crossval(X, Y, Tee, nK, K)
  
  # generate basis, add intercept
  X_poly <- cbind(1,polybasis(X,K_CV))
  
  # define M_P (I-P(P'P)^{-1}P)
  M_P <- diag(n) - (X_poly %*% solve((t(X_poly) %*% X_poly)) %*% t(X_poly))
  
  # estimate theta(K)
  theta <- (t(Tee) %*% M_P %*% Y)/(t(Tee) %*% M_P %*% Tee)
  
  # sigma (for variance estimate)
  sigma <- diag( as.numeric((M_P %*% (Y - Tee*as.numeric(theta))))^2 )
  
  # standard error
  bread <- solve((t(Tee) %*% M_P %*% Tee))
  se <- sqrt( bread %*% (t(Tee)%*%M_P%*%sigma%*%M_P%*%Tee) %*% bread)
  
  # save to vectors
  thetas[m] <- theta
  SEs[m]    <- se
  Ks[m]     <- K_CV
}
proc.time() - ptm
# runtime 14 minutes

# prep data for plots to show results
df <- data.frame(K = Ks, theta = thetas, se = SEs, rep = 1:M) %>%
  mutate(ci_low = theta - 1.96*se, ci_high = theta + 1.96*se) %>%
  arrange(ci_low) %>% mutate(rep_sorted = 1:M)

# panel 1: histogram of optimal Ks
k_df <- df %>% group_by(K) %>% summarise(Count = n()) %>% mutate(K = as.character(K))
p1 <- ggplot(k_df, aes(x=K,y=Count)) + geom_bar(stat='identity') +
  theme_minimal() + ggtitle('Optimal K')
p1

# panel 2: distribution of theta
p2 <- ggplot(df, aes(x=theta)) + geom_histogram(bins=12) + theme_minimal() +
  geom_vline(xintercept=1,linetype='solid',color='red',size=1) +
  labs(title='Estimates of Theta',x='Theta',y='Count')
p2

# panel 3: distribution of SEs
p3 <- ggplot(df, aes(x=se)) + geom_histogram(bins=12) + theme_minimal() +
  labs(title='SEs on Estimates of Theta',x='Standard Error',y='Count')
p3

# panel 4: confidence intervals, sorted by lower point
p4 <- ggplot(df) +
  geom_line(aes(x = rep_sorted, y=ci_low)) +
  geom_line(aes(x = rep_sorted, y=ci_high)) +
  geom_hline(yintercept=1,linetype='solid',color='red',size=1) +
  theme_minimal() +
  labs(title='Confidence Intervals for Theta',x='Replication (Sorted)',y='Theta')
p4

# combine with multiplot
source('multiplot.R')
png('q3_4c_R.png')
multiplot(p1, p3, p2, p4, cols=2)
dev.off()

# values for latex
avg_K     <- mean(Ks)
median_K  <- median(Ks)
avg_theta <- mean(thetas)
avg_bias  <- mean(thetas)-1
samp_var  <- sd(thetas)^2
avg_vhat  <- mean(SEs^2)
coverage  <- 1 - mean(df$ci_low > 1 | df$ci_high < 1)

###############################################################################