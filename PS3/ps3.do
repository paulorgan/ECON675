********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS3
* Last Update: Oct 22, 2018
********************************************************************************
clear all
set more off
capture log close

cd "C:\Users\prorgan\Box\Classes\Econ 675\Problem Sets\PS3"
log using ps3.log, replace

********************************************************************************
*** Question 1: Non-linear Least Squares
********************************************************************************
* Q1.9 setup
* load data
use pisofirme, clear

* add variable indicating having outcome data
gen s = 1-dmissing

* add log variable for use in regression [log(S_incomepc + 1)]
gen log_SincpcP1 = log(S_incomepc + 1)

* Q1.9a - estimates and statistics
* logit regression, robust standard errors
logit s S_age S_HHpeople log_SincpcP1, vce(robust)

* output for LaTeX
outreg2 using q1_9a.tex, side stats(coef se tstat pval ci) ///
 noaster noparen nor2 noobs dec(3) replace

* Q1.9b - nonparametric bootstrap
logit s S_age S_HHpeople log_SincpcP1, vce(bootstrap, reps(999))

* output for LaTeX
outreg2 using q1_9b.tex, side stats(coef se tstat pval ci) ///
 noaster noparen nor2 noobs dec(3) replace

* Q1.9c - propensity scores
* logit regression, robust standard errors
logit s S_age S_HHpeople log_SincpcP1, vce(robust)

* predict propensity score
predict p

* plot histogram, overlay kernel density
twoway histogram p || kdensity p, k(gaussian) || ///
 kdensity p, k(epanechnikov) || kdensity p, k(triangle) ///
 leg(lab(1 "Propensity Score") lab(2 "Gaussian") ///
	 lab(3 "Epanechnikov") lab(4 "Triangle"))
	 
* save
graph export q1_9c_S.png, replace

********************************************************************************
*** Question 2: Semiparametric GMM with Missing Data
********************************************************************************
* Q2.2b - feasible estimator
* since we assume F is logistic, we can drop the g_0 term
local vars = "dpisofirme S_age S_HHpeople log_SincpcP1"
gmm (danemia-invlogit({xb:`vars'})), instruments(`vars') vce(boot)

* output for LaTeX
outreg2 using q2_2b.tex, side stats(coef se tstat pval ci) ///
 noaster noparen nor2 noobs dec(3) replace

* Q2.3c - feasible estimator





*** other stuff from first attempts below, in case useful
* this runs and is close to my R results, but is not the moment condition we want
* but may actually be right? the g function is based on 

* need to somehow include the g_0 function in there, and account for conditional on s=1
local vars = "dpisofirme S_age S_HHpeople log_SincpcP1"
local f = logit({xb:`vars'})
local F = invlogit({xb:`vars'})
gmm ( (danemia - `F')*(`f'/(`F'*(1-`F'))) ), instruments(`vars') vce(bootstrap)

* doesn't work
* maybe define a program for g, then build it into the gmm moment condition?

* try writing it out:
local vars = "dpisofirme S_age S_HHpeople log_SincpcP1"
gmm ( (danemia - invlogit({xb:`vars'}))*(logit({xb:})/(invlogit({xb:})*(1-invlogit({xb:})))) ), instruments(`vars')

* also doesn't work

* might be able to use teffects?

********************************************************************************
*** Question 3: When Bootstrap Fails
********************************************************************************
* Q3.1 - nonparametric bootstrap
clear all

* generate sample
set seed 123
set obs 1000
gen X = runiform()

* save actual max
sum X
local maxX=r(max)

* run nonparametric bootstrap of max
bootstrap stat=r(max), reps(599) saving(nonpar_results, replace): summarize X

* load results
use nonpar_results, clear

* generate statistic
gen nonpar_stat = 1000*(`maxX'-stat)

* plot
hist nonpar_stat, ///
 plot(function exponential = 1-exponential(1,x), range(0 5) color(red))
graph export q3_1_S.png, replace

********************************************************************************
* Q3.2 - parametric bootstrap
clear all

tempname memhold
tempfile para_results

* generate sample
set seed 123
set obs 1000
gen X = runiform()

* save actual max
sum X
local maxX=r(max)

* parametric bootstrap
postfile `memhold' max using `para_results'
forvalues i = 1/599{
	capture drop sample
	gen sample = runiform(0,`maxX')
	sum sample
	post `memhold' (r(max))
}
postclose `memhold'

* load results
use `para_results', clear

* generate statistic
gen para_stat = 1000*(`maxX'-max)

* plot
hist para_stat, ///
 plot(function exponential = 1-exponential(1,x), range(0 5) color(red))
graph export q3_2_S.png, replace

********************************************************************************
log close
********************************************************************************
