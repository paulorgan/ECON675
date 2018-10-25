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

********************************************************************************
* Q2.1: Difference-in-Means

* experimental data
qui reg re78 t if (treat != 2), vce(hc2)

* PSID data
qui reg re78 t if (treat != 0), vce(hc2)

********************************************************************************
* Q2.2: Linear Least-Squares

* experimental data, three ways
qui reg re78 t `a' if (treat != 2), vce(r)
qui reg re78 t `b' if (treat != 2), vce(r)
qui reg re78 t `c' if (treat != 2), vce(r)

* PSID data
qui reg re78 t `a' if (treat != 0), vce(r)
qui reg re78 t `b' if (treat != 0), vce(r)
qui reg re78 t `c' if (treat != 0), vce(r)

********************************************************************************
* Q2.3: Regression Imputation

* ate, experimental data, three ways
qui teffects ra (re78 `a') (t) if (treat != 2), ate
qui teffects ra (re78 `b') (t) if (treat != 2), ate
qui teffects ra (re78 `c') (t) if (treat != 2), ate

* att, experimental data, three ways
qui teffects ra (re78 `a') (t) if (treat != 2), atet
qui teffects ra (re78 `b') (t) if (treat != 2), atet
qui teffects ra (re78 `c') (t) if (treat != 2), atet

* ate, PSID data, three ways
qui teffects ra (re78 `a') (t) if (treat != 0), ate
qui teffects ra (re78 `b') (t) if (treat != 0), ate
qui teffects ra (re78 `c') (t) if (treat != 0), ate

* att, PSID data, three ways
qui teffects ra (re78 `a') (t) if (treat != 0), atet
qui teffects ra (re78 `b') (t) if (treat != 0), atet
qui teffects ra (re78 `c') (t) if (treat != 0), atet

********************************************************************************
* Q2.4: Inverse Probability Weighting

* ate, experimental data, three ways
qui teffects ipw (re78) (t `a', probit) if (treat != 2), ate iterate(50)
qui teffects ipw (re78) (t `b', probit) if (treat != 2), ate iterate(50)
qui teffects ipw (re78) (t `c', probit) if (treat != 2), ate iterate(50)

* att, experimental data, three ways
qui teffects ipw (re78) (t `a', probit) if (treat != 2), atet iterate(50)
qui teffects ipw (re78) (t `b', probit) if (treat != 2), atet iterate(50)
qui teffects ipw (re78) (t `c', probit) if (treat != 2), atet iterate(50)

* predicted propensity scores for PSID data are too close to 0 or 1
* so we need to predict, then only use interior data
* note new if condition (keep if treated, or PSID and interior prop score)

* three ways, PSID data, ate and att
qui probit t `a' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipw (re78) (t `a', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
 
qui teffects ipw (re78) (t `a', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
 
qui probit t `b' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipw (re78) (t `b', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
 
qui teffects ipw (re78) (t `b', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
 
qui probit t `c' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipw (re78) (t `c', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)

qui teffects ipw (re78) (t `c', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)

********************************************************************************
* Q2.5: Doubly Robust

* ate, experimental data, three ways
qui teffects ipwra (re78 `a') (t `a', probit) if (treat != 2), ate iterate(50)
qui teffects ipwra (re78 `b') (t `b', probit) if (treat != 2), ate iterate(50)
qui teffects ipwra (re78 `c') (t `c', probit) if (treat != 2), ate iterate(50)

* att, experimental data, three ways
qui teffects ipwra (re78 `a') (t `a', probit) if (treat != 2), atet iterate(50)
qui teffects ipwra (re78 `b') (t `b', probit) if (treat != 2), atet iterate(50)
qui teffects ipwra (re78 `c') (t `c', probit) if (treat != 2), atet iterate(50)

* predicted propensity scores for PSID data are too close to 0 or 1
* so we need to predict, then only use interior data

* three ways, PSID data, ate and att
qui probit t `a' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipwra (re78 `a') (t `a', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
 
qui teffects ipwra (re78 `a') (t `a', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
 
qui probit t `b' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipwra (re78 `b') (t `b', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)
 
qui teffects ipwra (re78 `b') (t `b', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)
 
qui probit t `c' if (treat != 0)
capture: drop prop
predict prop
qui teffects ipwra (re78 `c') (t `c', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate iterate(50)

qui teffects ipwra (re78 `c') (t `c', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet iterate(50)

********************************************************************************
* Q2.6: Nearest Neighbor Matching

* ate, experimental data, three ways
qui teffects nnmatch (re78 `a') (t) if (treat != 2), ate nneighbor(1) metric(maha)
qui teffects nnmatch (re78 `b') (t) if (treat != 2), ate nneighbor(1) metric(maha)
qui teffects nnmatch (re78 `c') (t) if (treat != 2), ate nneighbor(1) metric(maha)

* att, experimental data, three ways
qui teffects nnmatch (re78 `a') (t) if (treat != 2), atet nneighbor(1) metric(maha)
qui teffects nnmatch (re78 `b') (t) if (treat != 2), atet nneighbor(1) metric(maha)
qui teffects nnmatch (re78 `c') (t) if (treat != 2), atet nneighbor(1) metric(maha)

* ate, PSID data, three ways
qui teffects nnmatch (re78 `a') (t) if (treat != 0), ate nneighbor(1) metric(maha)
qui teffects nnmatch (re78 `b') (t) if (treat != 0), ate nneighbor(1) metric(maha)
qui teffects nnmatch (re78 `c') (t) if (treat != 0), ate nneighbor(1) metric(maha)

* att, PSID data, three ways
qui teffects nnmatch (re78 `a') (t) if (treat != 0), atet nneighbor(1) metric(maha)
qui teffects nnmatch (re78 `b') (t) if (treat != 0), atet nneighbor(1) metric(maha)
qui teffects nnmatch (re78 `c') (t) if (treat != 0), atet nneighbor(1) metric(maha)

********************************************************************************
* Q2.7: Propensity Score Matching

* ate, experimental data, three ways
qui teffects psmatch (re78) (t `a', probit) if (treat != 2), ate
qui teffects psmatch (re78) (t `b', probit) if (treat != 2), ate
qui teffects psmatch (re78) (t `c', probit) if (treat != 2), ate

* att, experimental data, three ways
qui teffects psmatch (re78) (t `a', probit) if (treat != 2), atet
qui teffects psmatch (re78) (t `b', probit) if (treat != 2), atet
qui teffects psmatch (re78) (t `c', probit) if (treat != 2), atet

* three ways, PSID data, ate and att
qui probit t `a' if (treat != 0)
capture: drop prop
predict prop
qui teffects psmatch (re78) (t `a', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate
 
qui teffects psmatch (re78) (t `a', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet
 
qui probit t `b' if (treat != 0)
capture: drop prop
predict prop
qui teffects psmatch (re78) (t `b', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate
 
qui teffects psmatch (re78) (t `b', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet
 
qui probit t `c' if (treat != 0)
capture: drop prop
predict prop
qui teffects psmatch (re78) (t `c', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), ate

qui teffects psmatch (re78) (t `c', probit) ///
 if (treat==1 | (treat==2 & prop >= .0001 & prop <= .9999) ), atet

********************************************************************************
*** Question 3: Post-model Selection Inference
********************************************************************************
* Q3 setup

********************************************************************************
* Q3.1: Summary Statistics and Kernel Density Plots

********************************************************************************
* Q3.2: Empirical Coverage Rates

********************************************************************************
log close
********************************************************************************
