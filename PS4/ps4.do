********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS4
* Last Update: Oct 25, 2018
********************************************************************************
clear all
set more off
capture log close

cd "C:\Users\prorgan\Box\Classes\Econ 675\Problem Sets\PS4"
log using ps4.log, replace

********************************************************************************
*** Question 2: Estimating Average Treatment Effects
********************************************************************************
* Q2 setup

* load data
import delimited "LaLonde_all.csv"

* add treatment variable in expanded dataset
gen t = 0
replace t = 1 if treat == 1

* add other variables we will use in the analysis
gen logre74 = log(re74 + 1)
gen logre75 = log(re75 + 1)
gen age2 = age^2
gen educ2 = educ^2
gen age3 = age^3
gen bXu74 = black * u74
gen eXre74 = educ * logre74

* define vars for model specifications
local a = "age educ black hisp married nodegr logre74 logre75"
local b = "age educ black hisp married nodegr logre74 logre75 age2 educ2 u74 u75"
local c = "age educ black hisp married nodegr logre74 logre75 age2 educ2 u74 u75 age3 bXu74 eXre74"

* empty tables to fill with results (1+6*3)
matrix ate = J(19,8,.)
matrix att = J(19,8,.)

********************************************************************************
* Q2.1: Difference-in-Means (row 1)

* experimental data
qui reg re78 t if (treat != 2), vce(hc2)

* pull out results, and assign them to our empty table for ATE
matrix out = r(table)
matrix ate[1,1] = out["b","t"]
matrix ate[1,2] = out["se","t"]
matrix ate[1,3] = out["ll","t"]
matrix ate[1,4] = out["ul","t"]

* PSID data
qui reg re78 t if (treat != 0), vce(hc2)

* pull out results, and assign them to our empty table for ATE
matrix out = r(table)
matrix ate[1,5] = out["b","t"]
matrix ate[1,6] = out["se","t"]
matrix ate[1,7] = out["ll","t"]
matrix ate[1,8] = out["ul","t"]

* same for ATT
forvalues j = 1/8{
	matrix att[1,`j'] = ate[1,`j']
}

********************************************************************************
* Q2.2: Linear Least-Squares (rows 2-4)

* experimental data, three ways
qui reg re78 t `a' if (treat != 2), vce(r)
matrix out = r(table)
matrix ate[2,1] = out["b","t"]
matrix ate[2,2] = out["se","t"]
matrix ate[2,3] = out["ll","t"]
matrix ate[2,4] = out["ul","t"]

qui reg re78 t `b' if (treat != 2), vce(r)
matrix out = r(table)
matrix ate[3,1] = out["b","t"]
matrix ate[3,2] = out["se","t"]
matrix ate[3,3] = out["ll","t"]
matrix ate[3,4] = out["ul","t"]

qui reg re78 t `c' if (treat != 2), vce(r)
matrix out = r(table)
matrix ate[4,1] = out["b","t"]
matrix ate[4,2] = out["se","t"]
matrix ate[4,3] = out["ll","t"]
matrix ate[4,4] = out["ul","t"]

* PSID data
qui reg re78 t `a' if (treat != 0), vce(r)
matrix out = r(table)
matrix ate[2,5] = out["b","t"]
matrix ate[2,6] = out["se","t"]
matrix ate[2,7] = out["ll","t"]
matrix ate[2,8] = out["ul","t"]

qui reg re78 t `b' if (treat != 0), vce(r)
matrix out = r(table)
matrix ate[3,5] = out["b","t"]
matrix ate[3,6] = out["se","t"]
matrix ate[3,7] = out["ll","t"]
matrix ate[3,8] = out["ul","t"]

qui reg re78 t `c' if (treat != 0), vce(r)
matrix out = r(table)
matrix ate[4,5] = out["b","t"]
matrix ate[4,6] = out["se","t"]
matrix ate[4,7] = out["ll","t"]
matrix ate[4,8] = out["ul","t"]

* same for ATT
forvalues i = 2/4{
	forvalues j = 1/8{
		matrix att[`i',`j'] = ate[`i',`j']
	}
}

********************************************************************************
* Q2.3: Regression Imputation (rows 5-7)

* ate, experimental data, three ways
qui teffects ra (re78 `a') (t) if (treat != 2), ate
matrix out = r(table)
matrix ate[5,1] = out["b",1]
matrix ate[5,2] = out["se",1]
matrix ate[5,3] = out["ll",1]
matrix ate[5,4] = out["ul",1]

qui teffects ra (re78 `b') (t) if (treat != 2), ate
matrix out = r(table)
matrix ate[6,1] = out["b",1]
matrix ate[6,2] = out["se",1]
matrix ate[6,3] = out["ll",1]
matrix ate[6,4] = out["ul",1]

qui teffects ra (re78 `c') (t) if (treat != 2), ate
matrix out = r(table)
matrix ate[7,1] = out["b",1]
matrix ate[7,2] = out["se",1]
matrix ate[7,3] = out["ll",1]
matrix ate[7,4] = out["ul",1]

* att, experimental data, three ways
qui teffects ra (re78 `a') (t) if (treat != 2), atet
matrix out = r(table)
matrix att[5,1] = out["b",1]
matrix att[5,2] = out["se",1]
matrix att[5,3] = out["ll",1]
matrix att[5,4] = out["ul",1]

qui teffects ra (re78 `b') (t) if (treat != 2), atet
matrix out = r(table)
matrix att[6,1] = out["b",1]
matrix att[6,2] = out["se",1]
matrix att[6,3] = out["ll",1]
matrix att[6,4] = out["ul",1]

qui teffects ra (re78 `c') (t) if (treat != 2), atet
matrix out = r(table)
matrix att[7,1] = out["b",1]
matrix att[7,2] = out["se",1]
matrix att[7,3] = out["ll",1]
matrix att[7,4] = out["ul",1]

* ate, PSID data, three ways
qui teffects ra (re78 `a') (t) if (treat != 0), ate
matrix out = r(table)
matrix ate[5,5] = out["b",1]
matrix ate[5,6] = out["se",1]
matrix ate[5,7] = out["ll",1]
matrix ate[5,8] = out["ul",1]

qui teffects ra (re78 `b') (t) if (treat != 0), ate
matrix out = r(table)
matrix ate[6,5] = out["b",1]
matrix ate[6,6] = out["se",1]
matrix ate[6,7] = out["ll",1]
matrix ate[6,8] = out["ul",1]

qui teffects ra (re78 `c') (t) if (treat != 0), ate
matrix out = r(table)
matrix ate[7,5] = out["b",1]
matrix ate[7,6] = out["se",1]
matrix ate[7,7] = out["ll",1]
matrix ate[7,8] = out["ul",1]

* att, PSID data, three ways
qui teffects ra (re78 `a') (t) if (treat != 0), atet
matrix out = r(table)
matrix att[5,5] = out["b",1]
matrix att[5,6] = out["se",1]
matrix att[5,7] = out["ll",1]
matrix att[5,8] = out["ul",1]

qui teffects ra (re78 `b') (t) if (treat != 0), atet
matrix out = r(table)
matrix att[6,5] = out["b",1]
matrix att[6,6] = out["se",1]
matrix att[6,7] = out["ll",1]
matrix att[6,8] = out["ul",1]

qui teffects ra (re78 `c') (t) if (treat != 0), atet
matrix out = r(table)
matrix att[7,5] = out["b",1]
matrix att[7,6] = out["se",1]
matrix att[7,7] = out["ll",1]
matrix att[7,8] = out["ul",1]

********************************************************************************
* Q2.4: Inverse Probability Weighting (rows 8-10)

* ate, experimental data, three ways
qui teffects ipw (re78) (t `a', probit) if (treat != 2), ate iterate(50)
matrix out = r(table)
matrix ate[8,1] = out["b",1]
matrix ate[8,2] = out["se",1]
matrix ate[8,3] = out["ll",1]
matrix ate[8,4] = out["ul",1]

qui teffects ipw (re78) (t `b', probit) if (treat != 2), ate iterate(50)
matrix out = r(table)
matrix ate[9,1] = out["b",1]
matrix ate[9,2] = out["se",1]
matrix ate[9,3] = out["ll",1]
matrix ate[9,4] = out["ul",1]

qui teffects ipw (re78) (t `c', probit) if (treat != 2), ate iterate(50)
matrix out = r(table)
matrix ate[10,1] = out["b",1]
matrix ate[10,2] = out["se",1]
matrix ate[10,3] = out["ll",1]
matrix ate[10,4] = out["ul",1]

* att, experimental data, three ways
qui teffects ipw (re78) (t `a', probit) if (treat != 2), atet iterate(50)
matrix out = r(table)
matrix att[8,1] = out["b",1]
matrix att[8,2] = out["se",1]
matrix att[8,3] = out["ll",1]
matrix att[8,4] = out["ul",1]

qui teffects ipw (re78) (t `b', probit) if (treat != 2), atet iterate(50)
matrix out = r(table)
matrix att[9,1] = out["b",1]
matrix att[9,2] = out["se",1]
matrix att[9,3] = out["ll",1]
matrix att[9,4] = out["ul",1]

qui teffects ipw (re78) (t `c', probit) if (treat != 2), atet iterate(50)
matrix out = r(table)
matrix att[10,1] = out["b",1]
matrix att[10,2] = out["se",1]
matrix att[10,3] = out["ll",1]
matrix att[10,4] = out["ul",1]

* predicted propensity scores for PSID data are too close to 0 or 1
* so we need to predict, then only use interior data
* note new if condition (keep if treated, or PSID and interior prop score)

* three ways, PSID data, ate and att
qui probit t `a' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipw (re78) (t `a') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
matrix out = r(table)
matrix ate[8,5] = out["b",1]
matrix ate[8,6] = out["se",1]
matrix ate[8,7] = out["ll",1]
matrix ate[8,8] = out["ul",1]
 
qui teffects ipw (re78) (t `a') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
matrix out = r(table)
matrix att[8,5] = out["b",1]
matrix att[8,6] = out["se",1]
matrix att[8,7] = out["ll",1]
matrix att[8,8] = out["ul",1]
 
qui probit t `b' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipw (re78) (t `b') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
matrix out = r(table)
matrix ate[9,5] = out["b",1]
matrix ate[9,6] = out["se",1]
matrix ate[9,7] = out["ll",1]
matrix ate[9,8] = out["ul",1]
 
qui teffects ipw (re78) (t `b') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
matrix out = r(table)
matrix att[9,5] = out["b",1]
matrix att[9,6] = out["se",1]
matrix att[9,7] = out["ll",1]
matrix att[9,8] = out["ul",1]
 
qui probit t `c' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipw (re78) (t `c') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
matrix out = r(table)
matrix ate[10,5] = out["b",1]
matrix ate[10,6] = out["se",1]
matrix ate[10,7] = out["ll",1]
matrix ate[10,8] = out["ul",1]

qui teffects ipw (re78) (t `c') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
matrix out = r(table)
matrix att[10,5] = out["b",1]
matrix att[10,6] = out["se",1]
matrix att[10,7] = out["ll",1]
matrix att[10,8] = out["ul",1]

********************************************************************************
* Q2.5: Doubly Robust (rows 11-13)

* ate, experimental data, three ways
qui teffects ipwra (re78 `a') (t `a', probit) if (treat != 2), ate iterate(50)
matrix out = r(table)
matrix ate[11,1] = out["b",1]
matrix ate[11,2] = out["se",1]
matrix ate[11,3] = out["ll",1]
matrix ate[11,4] = out["ul",1]

qui teffects ipwra (re78 `b') (t `b', probit) if (treat != 2), ate iterate(50)
matrix out = r(table)
matrix ate[12,1] = out["b",1]
matrix ate[12,2] = out["se",1]
matrix ate[12,3] = out["ll",1]
matrix ate[12,4] = out["ul",1]

qui teffects ipwra (re78 `c') (t `c', probit) if (treat != 2), ate iterate(50)
matrix out = r(table)
matrix ate[13,1] = out["b",1]
matrix ate[13,2] = out["se",1]
matrix ate[13,3] = out["ll",1]
matrix ate[13,4] = out["ul",1]

* att, experimental data, three ways
qui teffects ipwra (re78 `a') (t `a', probit) if (treat != 2), atet iterate(50)
matrix out = r(table)
matrix att[11,1] = out["b",1]
matrix att[11,2] = out["se",1]
matrix att[11,3] = out["ll",1]
matrix att[11,4] = out["ul",1]

qui teffects ipwra (re78 `b') (t `b', probit) if (treat != 2), atet iterate(50)
matrix out = r(table)
matrix att[12,1] = out["b",1]
matrix att[12,2] = out["se",1]
matrix att[12,3] = out["ll",1]
matrix att[12,4] = out["ul",1]

qui teffects ipwra (re78 `c') (t `c', probit) if (treat != 2), atet iterate(50)
matrix out = r(table)
matrix att[13,1] = out["b",1]
matrix att[13,2] = out["se",1]
matrix att[13,3] = out["ll",1]
matrix att[13,4] = out["ul",1]

* predicted propensity scores for PSID data are too close to 0 or 1
* so we need to predict, then only use interior data

* three ways, PSID data, ate and att
qui probit t `a' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipwra (re78 `a') (t `a') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
matrix out = r(table)
matrix ate[11,5] = out["b",1]
matrix ate[11,6] = out["se",1]
matrix ate[11,7] = out["ll",1]
matrix ate[11,8] = out["ul",1]
 
qui teffects ipwra (re78 `a') (t `a') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
matrix out = r(table)
matrix att[11,5] = out["b",1]
matrix att[11,6] = out["se",1]
matrix att[11,7] = out["ll",1]
matrix att[11,8] = out["ul",1]
 
qui probit t `b' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipwra (re78 `b') (t `b') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
matrix out = r(table)
matrix ate[12,5] = out["b",1]
matrix ate[12,6] = out["se",1]
matrix ate[12,7] = out["ll",1]
matrix ate[12,8] = out["ul",1]
 
qui teffects ipwra (re78 `b') (t `b') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
matrix out = r(table)
matrix att[12,5] = out["b",1]
matrix att[12,6] = out["se",1]
matrix att[12,7] = out["ll",1]
matrix att[12,8] = out["ul",1]
 
qui probit t `c' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipwra (re78 `c') (t `c') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
matrix out = r(table)
matrix ate[13,5] = out["b",1]
matrix ate[13,6] = out["se",1]
matrix ate[13,7] = out["ll",1]
matrix ate[13,8] = out["ul",1]

qui teffects ipwra (re78 `c') (t `c') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
matrix out = r(table)
matrix att[13,5] = out["b",1]
matrix att[13,6] = out["se",1]
matrix att[13,7] = out["ll",1]
matrix att[13,8] = out["ul",1]

********************************************************************************
* Q2.6: Nearest Neighbor Matching (rows 14-16)

* ate, experimental data, three ways
qui teffects nnmatch (re78 `a') (t) if (treat != 2), ate nneighbor(1) metric(maha)
matrix out = r(table)
matrix ate[14,1] = out["b",1]
matrix ate[14,2] = out["se",1]
matrix ate[14,3] = out["ll",1]
matrix ate[14,4] = out["ul",1]

qui teffects nnmatch (re78 `b') (t) if (treat != 2), ate nneighbor(1) metric(maha)
matrix out = r(table)
matrix ate[15,1] = out["b",1]
matrix ate[15,2] = out["se",1]
matrix ate[15,3] = out["ll",1]
matrix ate[15,4] = out["ul",1]

qui teffects nnmatch (re78 `c') (t) if (treat != 2), ate nneighbor(1) metric(maha)
matrix out = r(table)
matrix ate[16,1] = out["b",1]
matrix ate[16,2] = out["se",1]
matrix ate[16,3] = out["ll",1]
matrix ate[16,4] = out["ul",1]

* att, experimental data, three ways
qui teffects nnmatch (re78 `a') (t) if (treat != 2), atet nneighbor(1) metric(maha)
matrix out = r(table)
matrix att[14,1] = out["b",1]
matrix att[14,2] = out["se",1]
matrix att[14,3] = out["ll",1]
matrix att[14,4] = out["ul",1]

qui teffects nnmatch (re78 `b') (t) if (treat != 2), atet nneighbor(1) metric(maha)
matrix out = r(table)
matrix att[15,1] = out["b",1]
matrix att[15,2] = out["se",1]
matrix att[15,3] = out["ll",1]
matrix att[15,4] = out["ul",1]

qui teffects nnmatch (re78 `c') (t) if (treat != 2), atet nneighbor(1) metric(maha)
matrix out = r(table)
matrix att[16,1] = out["b",1]
matrix att[16,2] = out["se",1]
matrix att[16,3] = out["ll",1]
matrix att[16,4] = out["ul",1]

* ate, PSID data, three ways
qui teffects nnmatch (re78 `a') (t) if (treat != 0), ate nneighbor(1) metric(maha)
matrix out = r(table)
matrix ate[14,5] = out["b",1]
matrix ate[14,6] = out["se",1]
matrix ate[14,7] = out["ll",1]
matrix ate[14,8] = out["ul",1]

qui teffects nnmatch (re78 `b') (t) if (treat != 0), ate nneighbor(1) metric(maha)
matrix out = r(table)
matrix ate[15,5] = out["b",1]
matrix ate[15,6] = out["se",1]
matrix ate[15,7] = out["ll",1]
matrix ate[15,8] = out["ul",1]

qui teffects nnmatch (re78 `c') (t) if (treat != 0), ate nneighbor(1) metric(maha)
matrix out = r(table)
matrix ate[16,5] = out["b",1]
matrix ate[16,6] = out["se",1]
matrix ate[16,7] = out["ll",1]
matrix ate[16,8] = out["ul",1]

* att, PSID data, three ways
qui teffects nnmatch (re78 `a') (t) if (treat != 0), atet nneighbor(1) metric(maha)
matrix out = r(table)
matrix att[14,5] = out["b",1]
matrix att[14,6] = out["se",1]
matrix att[14,7] = out["ll",1]
matrix att[14,8] = out["ul",1]

qui teffects nnmatch (re78 `b') (t) if (treat != 0), atet nneighbor(1) metric(maha)
matrix out = r(table)
matrix att[15,5] = out["b",1]
matrix att[15,6] = out["se",1]
matrix att[15,7] = out["ll",1]
matrix att[15,8] = out["ul",1]

qui teffects nnmatch (re78 `c') (t) if (treat != 0), atet nneighbor(1) metric(maha)
matrix out = r(table)
matrix att[16,5] = out["b",1]
matrix att[16,6] = out["se",1]
matrix att[16,7] = out["ll",1]
matrix att[16,8] = out["ul",1]

********************************************************************************
* Q2.7: Propensity Score Matching (rows 17-19)

* ate, experimental data, three ways
qui teffects psmatch (re78) (t `a', probit) if (treat != 2), ate
matrix out = r(table)
matrix ate[17,1] = out["b",1]
matrix ate[17,2] = out["se",1]
matrix ate[17,3] = out["ll",1]
matrix ate[17,4] = out["ul",1]

qui teffects psmatch (re78) (t `b', probit) if (treat != 2), ate
matrix out = r(table)
matrix ate[18,1] = out["b",1]
matrix ate[18,2] = out["se",1]
matrix ate[18,3] = out["ll",1]
matrix ate[18,4] = out["ul",1]

qui teffects psmatch (re78) (t `c', probit) if (treat != 2), ate
matrix out = r(table)
matrix ate[19,1] = out["b",1]
matrix ate[19,2] = out["se",1]
matrix ate[19,3] = out["ll",1]
matrix ate[19,4] = out["ul",1]

* att, experimental data, three ways
qui teffects psmatch (re78) (t `a', probit) if (treat != 2), atet
matrix out = r(table)
matrix att[17,1] = out["b",1]
matrix att[17,2] = out["se",1]
matrix att[17,3] = out["ll",1]
matrix att[17,4] = out["ul",1]

qui teffects psmatch (re78) (t `b', probit) if (treat != 2), atet
matrix out = r(table)
matrix att[18,1] = out["b",1]
matrix att[18,2] = out["se",1]
matrix att[18,3] = out["ll",1]
matrix att[18,4] = out["ul",1]

qui teffects psmatch (re78) (t `c', probit) if (treat != 2), atet
matrix out = r(table)
matrix att[19,1] = out["b",1]
matrix att[19,2] = out["se",1]
matrix att[19,3] = out["ll",1]
matrix att[19,4] = out["ul",1]

* three ways, PSID data, ate and att
qui probit t `a' if (treat != 0)
capture: drop prop
predict prop
qui teffects psmatch (re78) (t `a') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate
matrix out = r(table)
matrix ate[17,5] = out["b",1]
matrix ate[17,6] = out["se",1]
matrix ate[17,7] = out["ll",1]
matrix ate[17,8] = out["ul",1]
 
qui teffects psmatch (re78) (t `a') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet
matrix out = r(table)
matrix att[17,5] = out["b",1]
matrix att[17,6] = out["se",1]
matrix att[17,7] = out["ll",1]
matrix att[17,8] = out["ul",1]
 
qui probit t `b' if (treat != 0)
capture: drop prop
predict prop
qui teffects psmatch (re78) (t `b') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate
matrix out = r(table)
matrix ate[18,5] = out["b",1]
matrix ate[18,6] = out["se",1]
matrix ate[18,7] = out["ll",1]
matrix ate[18,8] = out["ul",1]
 
qui teffects psmatch (re78) (t `b') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet
matrix out = r(table)
matrix att[18,5] = out["b",1]
matrix att[18,6] = out["se",1]
matrix att[18,7] = out["ll",1]
matrix att[18,8] = out["ul",1]
 
qui probit t `c' if (treat != 0)
capture: drop prop
predict prop
qui teffects psmatch (re78) (t `c') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate
matrix out = r(table)
matrix ate[19,5] = out["b",1]
matrix ate[19,6] = out["se",1]
matrix ate[19,7] = out["ll",1]
matrix ate[19,8] = out["ul",1]

qui teffects psmatch (re78) (t `c') ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet
matrix out = r(table)
matrix att[19,5] = out["b",1]
matrix att[19,6] = out["se",1]
matrix att[19,7] = out["ll",1]
matrix att[19,8] = out["ul",1]

********************************************************************************
* Q2: Output table to Excel, then convert that to LaTeX

* write to c2 so we can add row, heading columns in Excel
putexcel set ate_s, modify
putexcel b3 = matrix(ate)
putexcel close

putexcel set att_s, modify
putexcel b3 = matrix(att)
putexcel close

********************************************************************************
*** Question 3: Post-model Selection Inference
********************************************************************************
* Q3 setup
clear all

********************************************************************************
* Q3.1: Summary Statistics and Kernel Density Plots

********************************************************************************
* Q3.2: Empirical Coverage Rates

********************************************************************************
log close
********************************************************************************
