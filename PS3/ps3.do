********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS3
* Last Update: Oct 26, 2018
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
* Q2.2b(MCAR)- feasible estimator

* gmm, four moment conditions
local vars = "dpisofirme S_age S_HHpeople log_SincpcP1"
gmm ((danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*dpisofirme) ///
((danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*S_age) ///
((danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*S_HHpeople) ///
((danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*log_SincpcP1), ///
instruments(`vars') winitial(identity) vce(boot)
 
* output for LaTeX
 mata: 
	coef = st_matrix("e(b)")'
	se = st_matrix("e(se)")'
	
	tstat = coef:/se
		
	CI_low = coef - 1.96:*se
	CI_high = coef + 1.96:*se
	 
	stats = round((coef,se,tstat,CI_low,CI_high),.001)
		
	st_matrix("stats",stats)
end
mat rownames stats = `vars'
mat colnames stats = coef se tstat CI_low CI_high
outtable using q2_2b, mat(stats) replace nobox

* Q2.3c (MAR)- feasible estimator
* we predicted p before, but did not use t, so do that now:
* logit regression, robust standard errors
logit s dpisofirme S_age S_HHpeople log_SincpcP1, vce(robust)

* predict propensity score
predict p_witht

* now run gmm adding in new term s/p
local vars = "dpisofirme S_age S_HHpeople log_SincpcP1"
gmm ((s/p_witht)*(danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*dpisofirme) ///
((s/p_witht)*(danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*S_age) ///
((s/p_witht)*(danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*S_HHpeople) ///
((s/p_witht)*(danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*log_SincpcP1), ///
instruments(`vars') winitial(identity) vce(boot)

* output for LaTeX
 mata: 
	coef = st_matrix("e(b)")'
	se = st_matrix("e(se)")'
	
	tstat = coef:/se
		
	CI_low = coef - 1.96:*se
	CI_high = coef + 1.96:*se
	 
	stats = round((coef,se,tstat,CI_low,CI_high),.001)
		
	st_matrix("stats",stats)
end
mat rownames stats = `vars'
mat colnames stats = coef se tstat CI_low CI_high
outtable using q2_3c, mat(stats) replace nobox

* Q2.3d (MAR)- feasible estimator, trimmed
* we predicted p before, and have s, so add that before the moment conditions
local vars = "dpisofirme S_age S_HHpeople log_SincpcP1"
gmm ((s/p_witht)*(danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*dpisofirme) ///
((s/p_witht)*(danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*S_age) ///
((s/p_witht)*(danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*S_HHpeople) ///
((s/p_witht)*(danemia - invlogit((dpisofirme*{theta}+S_age*{gamma1}+S_HHpeople*{gamma2}+log_SincpcP1*{gamma3})))*log_SincpcP1) ///
if p_witht >= 0.1, instruments(`vars') winitial(identity) vce(boot)

* output for LaTeX
 mata: 
	coef = st_matrix("e(b)")'
	se = st_matrix("e(se)")'
	
	tstat = coef:/se
		
	CI_low = coef - 1.96:*se
	CI_high = coef + 1.96:*se
	 
	stats = round((coef,se,tstat,CI_low,CI_high),.001)
		
	st_matrix("stats",stats)
end
mat rownames stats = `vars'
mat colnames stats = coef se tstat CI_low CI_high
outtable using q2_3d, mat(stats) replace nobox

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
