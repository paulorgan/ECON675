*******************************************************************
* ECON 675, Assignment 5
* Fall 2015
* University of Michigan
* Latest update: Oct 27, 2015
*******************************************************************

*******************************************************************
* Question 2
*******************************************************************
clear all
set more off
cap log close

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
    reps(5000) seed(123) nodots: ///
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

mat2txt, matrix(Results) saving(result1.txt) format(%9.4f) replace



*******************************************************************
* Question 3
*******************************************************************
clear all
set more off
cap log close
use "Angrist_Krueger.dta"

*************************************************************************
* The following replicates Columns (5)-(8), Table V 
* in Angrist and Krueger (1991 QJE)
*************************************************************************

*** Column 5, Table V, Angrist and Krueger (1991 QJE)
reg l_w_wage educ non_white married SMSA i.region i.YoB_ld

*** Column 6, Table V, Angrist and Krueger (1991 QJE)
ivregress 2sls l_w_wage non_white married SMSA i.region i.YoB_ld ///
    (educ = i.YoB_ld##i.QoB)
estat firststage

*** Column 7, Table V, Angrist and Krueger (1991 QJE)
reg  l_w_wage educ non_white married SMSA age_q age_sq i.region i.YoB_ld  

*** Column 8, Table V, Angrist and Krueger (1991 QJE)
ivregress 2sls l_w_wage non_white married SMSA age_q age_sq i.region i.YoB_ld ///
    (educ = i.YoB_ld##i.QoB)
estat firststage

*************************************************************************
* The following replicates Columns (1) and (2), Table 3 
* in Bound et al. (1995)
*************************************************************************
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


permute QoB TSLS_1_b = r(beta), reps(500) seed(123) saving(premute1, replace): ///
    IV_quick QoB, model(1)
	
permute QoB TSLS_2_b = _b[educ], reps(500) seed(123) saving(premute2, replace): ///
    ivregress 2sls l_w_wage non_white married SMSA age_q age_sq i.region i.YoB_ld ///
    (educ = i.YoB_ld##i.QoB)

clear all
use "premute1.dta"
sum TSLS_1_b

clear all
use "premute2.dta"
sum TSLS_2_b




