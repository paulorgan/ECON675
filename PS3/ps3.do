********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS3
* Last Update: Oct 19, 2018
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

* Q1.9a - estimates and statistics

* Q1.9b - nonparametric bootstrap

* Q1.9c - propensity scores

********************************************************************************
*** Question 2: Semiparametric GMM with Missing Data
********************************************************************************
* Q2.2b - feasible estimator

* Q2.3c - feasible estimator

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
