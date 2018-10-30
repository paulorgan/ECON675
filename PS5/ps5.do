********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS5
* Last Update: Oct 30, 2018
********************************************************************************
clear all
set more off
capture log close

cd "C:\Users\prorgan\Box\Classes\Econ 675\Problem Sets\PS5"
log using ps5.log, replace

********************************************************************************
*** Question 2: Weak Instrument Simulations
********************************************************************************
* Q2 setup

********************************************************************************
* Q2: OLS


********************************************************************************
* Q2: 2SLS


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
local iv_short = "non_white SMSA married"
local ols_long = "educ non_white SMSA married age_q age_sq"
local iv_long  = "non_white SMSA married age_q age_sq"

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

* somehow use the permute command (we used it in PS1 with this code:)
* Fisher permutation
* permute treat diffmean=(r(mu_2)-r(mu_1)), reps(999) nowarn: ttest earn78, by(treat)

local iv_short = "non_white SMSA married"
permute QoB beta = _b["educ"], reps(25) seed(22): ///
	ivregress2 2sls l_w_wage `iv_short' i.region i.YoB_ld ///
	(educ = i.YoB_ld##i.QoB)

* not sure this is best way to do it
* maybe store results along the way, then collapse?
	
********************************************************************************
log close
********************************************************************************
