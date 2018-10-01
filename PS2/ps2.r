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
require(ggplot2)   # plots
require(kedd)      # kernel estimation

set.seed(22)
select = dplyr::select
setwd('C:/Users/prorgan/Box/Classes/Econ 675/Problem Sets/PS2')

###############################################################################
## Question 1: Kernel Density Estimation
# Q1.3a

# sample size
n     <- 1000

# equally weight two distributions
comps <- sample(1:2,prob=c(.5,.5),size=n, replace=T)

# Normal density specs
mus <- c(-1.5, 1)
sds <- sqrt(c(1.5, 1))

# generate sample
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

# define DGP


# compute optimal bandwidth
h.amise(kernel = 'epanechnikov', deriv.order= 0)

# Q1.3b
# simulate 1000 times
M <- 1000

###############################################################################