********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS6
* Last Update: Nov 26, 2018
********************************************************************************
clear all
set more off
capture log close

cd "C:\Users\prorgan\Box\Classes\Econ 675\Problem Sets\PS6"
log using ps6.log, replace

use HeadStart.dta

* from the appendix
* counties are treated if povrate60 >0, not treated o/w
* outcome var is mort_related_post
* preintervention var is more_related_pre
* exog post var is mort_injury_post

* packages: rdrobust, rdlocrand, rddensity, rdmulti, rdpower
* https://sites.google.com/site/rdpackages/home

********************************************************************************
*** Question 2: The Effect of Head Start on Child Mortality
********************************************************************************
* Q2.1 (a): RD Plots

* using the rdrobust package from Matias
rdplot mort_related_pre povrate60, c(0) binselect(es) ///
    graph_options(title("Evenly-spaced binning, IMSE-optimal")) saving("s\1_1a")
graph export "s\1_1a.png", replace

rdplot mort_related_pre povrate60, c(0) binselect(esmv) ///
    graph_options(title("Evenly-spaced binning, variability-mimicking")) saving("s\1_1b")
graph export "s\1_1b.png", replace

rdplot mort_related_pre povrate60, c(0) binselect(qs) ///
    graph_options(title("Quantile-spaced binning, IMSE-optimal"))
graph export "s\1_1c.png", replace

rdplot mort_related_pre povrate60, c(0) binselect(qsmv) ///
    graph_options(title("Quantile-spaced binning, variability mimicking"))
graph export "s\1_1d.png", replace

********************************************************************************
* Q2.1 (b): Falsification Tests
* falsification of RD design using three methods

* (1) histogram plots
twoway (hist povrate60 if povrate60 <=0, freq bcolor(black)) ///
	   (hist povrate60 if povrate60 >0, freq bcolor(red))
graph export "s\1_2i.png", replace
	   
* (2) binomial tests using rdlocrand package - how to interpret?
rdwinselect povrate60, wmin(0.05) wstep(0.05) nwindows(100)

* (3) continuity-in-density tests using rddensity package - how to interpret?
rddensity povrate60, all
rddensity povrate60, all h(5 5)
rddensity povrate60, all h(1 1)

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
