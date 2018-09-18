********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS1
* Last Update: Sept 18, 2018
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
local V0 = r(Var)/(`N0'-1)

sum earn78 if treat==1
local N1 = r(N)
local mu1 = r(mean)
local sd1 = r(sd)
local V1 = r(Var)/(`N1'-1)

local tau = `mu1'-`mu0'
local v = sqrt(`V1'+`V0')
local T = `tau'/`v'
local pval = 2*normal(-abs(`T'))

local ci_low  = `tau' + invnormal(.025)*`v'
local ci_high = `tau' + invnormal(.975)*`v'

* rounding
local mu0 = round(`mu0', .01)
local mu1 = round(`mu1', .0001)
local sd0 = round(`sd0', .01)
local sd1 = round(`sd1', .0001)
local T = round(`T', .01)
local pval = round(`pval', .0001)
local ci_low = round(`ci_low', .01)
local ci_high = round(`ci_high', .01)

di "Control: `mu0' (`sd0') -- N=`N0'" 
di "Treatment: `mu1' (`sd1') -- N=`N1'" 
di "Test = `T' -- p-val = `pval'"
di "CI: [`ci_low', `ci_high']"

* canned version for comparison:
ttest earn78, by(treat) unequal

* 2) Fisher's Approach
* 2a) p-Value
* Fisher permutation
permute treat diffmean=(r(mu_2)-r(mu_1)), reps(999) nowarn: ttest earn78, by(treat)

* Kolgomorov-Smirnov
ksmirnov earn78, by(treat) exact

* 2b) confidence interval
* define imputed variables
gen y1_imp = earn78
replace y1_imp = earn78 + `tau' if treat==0

gen y0_imp = earn78
replace y0_imp = earn78 - `tau' if treat==1

* define program to bootstrap
program define SUTVAdiff, rclass
	sum y1_imp if treat==1
	local temp1 = r(mean)
	sum y0_imp if treat==0
	local temp0 = r(mean)
	return scalar SUTVAdiff = `temp1' - `temp0'
end

bootstrap diff = r(SUTVAdiff), reps(999): SUTVAdiff

* 3) Power Calculations
* 3a) graphing power function
local a = 0.05
local Z = invnormal(1-`a'/2)

** Plot power functions
twoway (function y = 1 - normal(x/`v'+`Z') + normal(x/`v'-`Z'), range(-2500 2500)), ///
	   yline(`a', lpattern(dash))
graph save power_stata.png, replace

* 3b) determining minimum sample size
mata:
	y = st_data(., "earn78"); t = st_data(., "treat")
	
	p = 2/3
	Z = 1.96
	tau = 1000
	
	V1 = quadvariance(select(y,t:==1))
	V0 = quadvariance(select(y,t:==0))
	
	N = (0::980)*5 + J(981,1,100)
	
	N1 = N*p
	N0 = N-N1
	
	std_errors = J(981,1,0)
	for(i = 1; i<=981; i++){
		std_errors[i] = sqrt(V1/N1[i]+V0/N0[i])
	}
	
	powers = J(981,1,0)
	for(i = 1; i<=981; i++){
		powers[i] = 1-normal(tau/std_errors[i]+Z)+normal(tau/std_errors[i]-Z)
	}
end

* display sample sizes and corresponding power to detect 1000 increase in earnings
mata: N, powers
********************************************************************************
