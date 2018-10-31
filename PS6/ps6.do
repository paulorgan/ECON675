********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS6
* Last Update: Oct 31, 2018
********************************************************************************
clear all
set more off
capture log close

cd "C:\Users\prorgan\Box\Classes\Econ 675\Problem Sets\PS6"
log using ps6.log, replace

********************************************************************************
*** Question 2: The Effect of Head Start on Child Mortality
********************************************************************************
* Q2.1 (a): RD Plots

* four plots: {evenly-spaced, quantile-spaced} X (IMSE-optimal,data var}

********************************************************************************
* Q2.1 (b): Falsification Tests
* falsification of RD design using three methods

* (1) histogram plots

* (2) binomial tests

* (3) continuity-in-density tests

********************************************************************************
* Q2.2 (a): constant treatment effect model

********************************************************************************
* Q2.2 (b): heterogeneous treatment effect model

********************************************************************************
* Q2.2 (c): local parametric model

********************************************************************************
* Q2.3 (a): MSE-optimal RD estimators

********************************************************************************
* Q2.3 (b): Robustness checks

********************************************************************************
* Q2.4: Local Randomization Methods
	
********************************************************************************
log close
********************************************************************************
