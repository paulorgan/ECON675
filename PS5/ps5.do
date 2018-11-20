********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS5
* Last Update: Nov 19, 2018
********************************************************************************
clear all
set more off
capture log close

cd "C:\Users\prorgan\Box\Classes\Econ 675\Problem Sets\PS5"
log using ps5.log, replace

********************************************************************************
*** Question 2: Weak Instrument Simulations
********************************************************************************

* from Yingjie
program define weak_IV, rclass
    syntax [, obs(integer 200) f_stat(real 10) ]
	drop _all
	
	set obs `obs'
	
	* DGP
	gen u = rnormal()
	gen v = 0.99 * u + sqrt(1-0.99^2) * rnormal()
	gen z = rnormal()
	
	local gamma_0 = sqrt((`f_stat' - 1) / `obs')
	gen x = `gamma_0' * z + v
	gen y = u
	
	* OLS
	qui reg y x, robust
	return scalar OLS_b   = _b[x]
	return scalar OLS_se  = _se[x]
	return scalar OLS_rej = abs(_b[x]/_se[x]) > 1.96
	
	* 2SLS
	qui ivregress 2sls y (x = z)
	return scalar TSLS_b   = _b[x]
	return scalar TSLS_se  = _se[x]
	return scalar TSLS_rej = abs(_b[x]/_se[x]) > 1.96
	qui reg x z
	return scalar TSLS_F   = e(F)
end

* simulation 1: F = 1 
simulate OLS_b=r(OLS_b) OLS_se=r(OLS_se) OLS_rej=r(OLS_rej) ///
    TSLS_b=r(TSLS_b) TSLS_se=r(TSLS_se) TSLS_rej=r(TSLS_rej) TSLS_F=r(TSLS_F), ///
    reps(5000) seed(22) nodots: ///
    weak_IV, f_stat(1)
	
local k = 1
matrix Results = J(7, 5, .)

qui sum OLS_b, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum OLS_se, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum OLS_rej, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_b, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_se, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_rej, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_F, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

mat2txt, matrix(Results) saving(q2_1.txt) format(%9.4f) replace

* now run for F=1.25, 10, and 100

* simulation 2: F = 1.25 
simulate OLS_b=r(OLS_b) OLS_se=r(OLS_se) OLS_rej=r(OLS_rej) ///
    TSLS_b=r(TSLS_b) TSLS_se=r(TSLS_se) TSLS_rej=r(TSLS_rej) TSLS_F=r(TSLS_F), ///
    reps(5000) seed(22) nodots: ///
    weak_IV, f_stat(1.25)
	
local k = 1
matrix Results = J(7, 5, .)

qui sum OLS_b, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum OLS_se, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum OLS_rej, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_b, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_se, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_rej, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_F, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

mat2txt, matrix(Results) saving(q2_2.txt) format(%9.4f) replace

* simulation 3: F = 10
simulate OLS_b=r(OLS_b) OLS_se=r(OLS_se) OLS_rej=r(OLS_rej) ///
    TSLS_b=r(TSLS_b) TSLS_se=r(TSLS_se) TSLS_rej=r(TSLS_rej) TSLS_F=r(TSLS_F), ///
    reps(5000) seed(22) nodots: ///
    weak_IV, f_stat(10)
	
local k = 1
matrix Results = J(7, 5, .)

qui sum OLS_b, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum OLS_se, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum OLS_rej, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_b, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_se, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_rej, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_F, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

mat2txt, matrix(Results) saving(q2_3.txt) format(%9.4f) replace

* simulation 4: F = 100 
simulate OLS_b=r(OLS_b) OLS_se=r(OLS_se) OLS_rej=r(OLS_rej) ///
    TSLS_b=r(TSLS_b) TSLS_se=r(TSLS_se) TSLS_rej=r(TSLS_rej) TSLS_F=r(TSLS_F), ///
    reps(5000) seed(22) nodots: ///
    weak_IV, f_stat(100)
	
local k = 1
matrix Results = J(7, 5, .)

qui sum OLS_b, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum OLS_se, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum OLS_rej, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_b, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_se, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_rej, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

qui sum TSLS_F, detail
matrix Results[`k',1] = r(mean)
matrix Results[`k',2] = r(sd)
matrix Results[`k',3] = r(p10)
matrix Results[`k',4] = r(p50)
matrix Results[`k',5] = r(p90)
local k = `k' + 1

mat2txt, matrix(Results) saving(q2_4.txt) format(%9.4f) replace

********************************************************************************
*** Question 3: Weak Instrument - Empirical Study
********************************************************************************
* Q3 setup
clear all

use Angrist_Krueger

********************************************************************************
* Q3.1: Angrist and Krueger (1991)

* variables for use in regressions
local ols_short = "educ non_white SMSA married"
local iv_short  = "non_white SMSA married"
local ols_long  = "educ non_white SMSA married age_q age_sq"
local iv_long   = "non_white SMSA married age_q age_sq"

* note going out of order from PS to match with A&K(1991)

* OLS 1 (these SEs match the table, so I don't think they use robust SEs)
reg l_w_wage `ols_short' i.region i.YoB_ld

* output for LaTeX
outreg2 using q3_1.tex, stats(coef se) keep(`ols_short') noaster dec(4) replace ///
	addtext(9 Year-of-birth dummies, Yes, 8 Region-of-residence dummies, Yes) ///
	ctitle(OLS 1) nor2 noobs

* 2SLS 1
ivregress2 2sls l_w_wage `iv_short' i.region i.YoB_ld ///
	(educ = i.YoB_ld##i.QoB)

* output for LaTeX
outreg2 using q3_1.tex, stats(coef se) keep(`ols_short') noaster dec(4) append ///
	addtext(9 Year-of-birth dummies, Yes, 8 Region-of-residence dummies, Yes) ///
	ctitle(2SLS 1) nor2 noobs
			
* OLS 2
reg l_w_wage `ols_long' i.region i.YoB_ld

* output for LaTeX
outreg2 using q3_1.tex, stats(coef se) keep(`ols_long') noaster dec(4) append ///
	addtext(9 Year-of-birth dummies, Yes, 8 Region-of-residence dummies, Yes) ///
	ctitle(OLS 2) nor2 noobs

* 2SLS 2
ivregress2 2sls l_w_wage `iv_long' i.region i.YoB_ld ///
    (educ = i.YoB_ld##i.QoB)
	
* output for LaTeX
outreg2 using q3_1.tex, stats(coef se) keep(`ols_long') noaster dec(4) append ///
	addtext(9 Year-of-birth dummies, Yes, 8 Region-of-residence dummies, Yes) ///
	ctitle(2SLS 2) nor2 noobs

********************************************************************************
* Q3.2: Bound, Jaeger, and Baker (1995)

* program defined in the problem set appendix for quicker IV estimation
capture program drop IV_quick
program define IV_quick, rclass
    syntax varlist(max=1) [, model(integer 1) ]
	local x "`varlist'"
	
	if (`model' == 1) {
	    capture drop educ_hat
		qui reg educ non_white married SMSA i.region i.YoB_ld i.YoB_ld##i.`x'
		predict educ_hat
		qui reg l_w_wage educ_hat non_white married SMSA i.region i.YoB_ld
		return scalar beta = _b[educ_hat]
	}
	if (`model' == 2) {
	    capture drop educ_hat
		qui reg educ non_white married SMSA age_q age_sq i.region i.YoB_ld i.YoB_ld##i.`x'
		predict educ_hat
		qui reg l_w_wage educ_hat non_white married SMSA age_q age_sq i.region i.YoB_ld
		return scalar beta = _b[educ_hat]
	}
end

* permute QoB 500 times for each model, save results
permute QoB TSLS_1_b = r(beta), reps(500) seed(22) saving(q3_model1, replace): ///
    IV_quick QoB, model(1)
	
permute QoB TSLS_2_b = r(beta), reps(500) seed(22) saving(q3_model2, replace): ///
    IV_quick QoB, model(2)

* summarize results necessary for table in LaTeX
clear all
use q3_model1
sum TSLS_1_b

clear all
use q3_model2
sum TSLS_2_b
	
********************************************************************************
log close
********************************************************************************
