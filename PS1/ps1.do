********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS1
* Last Update: Sept 12, 2018
********************************************************************************
clear all
set more off
capture log close

cd "C:\Users\prorgan\Box\Classes\Econ 675\Problem Sets\PS1"

log using ps1.log, replace

insheet using "LaLonde_1986.csv", clear

********************************************************************************
* Question 2: Implementing Least-Squares Estimators
* generate new variables
gen educ2 = educ^2
gen blackXearn74 = black * earn74

* run regression
reg earn78 treat black age educ educ2 earn74 blackXearn74 u74 u75

* output table to include in LaTeX
outreg2 using q2outreg.tex, side stats(coef se tstat pval ci) ///
 noaster noparen nor2 noobs dec(2) replace

********************************************************************************
* Question 3: Analysis of Experiments
* 1) Neyman's Approach
* manually
sum earn78 if treat==0
local N0 = r(N)
local mu0 = r(mean)
local sd0 = r(sd)
local V0 = r(Var)/r(N)

sum earn78 if treat==1
local N1 = r(N)
local mu1 = r(mean)
local sd1 = r(sd)
local V1 = r(Var)/r(N)

local tau = `mu1'-`mu0'
local v = sqrt(`V1'+`V0')
local T = `tau'/`v'
local pval = 2*normal(-abs(`T'))

local mu0 = round(`mu0', .01)
local mu1 = round(`mu1', .0001)
local sd0 = round(`sd0', .01)
local sd1 = round(`sd1', .0001)

di "Control: `mu0' (`sd0') -- N=`N0'" 
di "Treatment: `mu1' (`sd1') -- N=`N1'" 
di "Test = `T' -- p-val = `pval'"

* canned version for comparison:
ttest earn78, by(treat) unequal

* 2) Fisher's Approach
* Fisher permutation
permute treat diffmean=(r(mu_2)-r(mu_1)), reps(999) nowarn: ttest earn78, by(treat)

* Kolgomorov-Smirnov
ksmirnov earn78, by(treat) exact

* 3) Power Calculations


********************************************************************************
