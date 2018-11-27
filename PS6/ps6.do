********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS6
* Last Update: Nov 27, 2018
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

* rename variable for ease of use
rename povrate60 pov
rename mort_related_pre rel_pre
rename mort_related_post rel_post
rename mort_injury_post inj_post

********************************************************************************
*** Question 2: The Effect of Head Start on Child Mortality
********************************************************************************
* Q2.1.1: RD Plots

* using the rdrobust package from Matias
rdplot rel_pre pov, c(0) binselect(es) ///
    graph_options(title("Evenly-spaced binning, IMSE-optimal"))
graph export "s\1_1a.png", replace

rdplot rel_pre pov, c(0) binselect(esmv) ///
    graph_options(title("Evenly-spaced binning, variability-mimicking"))
graph export "s\1_1b.png", replace

rdplot rel_pre pov, c(0) binselect(qs) ///
    graph_options(title("Quantile-spaced binning, IMSE-optimal"))
graph export "s\1_1c.png", replace

rdplot rel_pre pov, c(0) binselect(qsmv) ///
    graph_options(title("Quantile-spaced binning, variability mimicking"))
graph export "s\1_1d.png", replace

********************************************************************************
* Q2.1.2: Falsification Tests
* falsification of RD design using three methods

gen t = pov > 0

* (1) histogram plots
twoway (hist pov if t, freq bcolor(black)) ///
	   (hist pov if !t, freq bcolor(red))
graph export "s\1_2i.png", replace
	   
* (2) binomial tests using rdlocrand package - how to interpret?
rdwinselect pov, wmin(0.05) wstep(0.05) nwindows(100)

* (3) continuity-in-density tests using rddensity package - how to interpret?
rddensity pov, all
rddensity pov, all h(5 5)
rddensity pov, all h(1 1)

********************************************************************************
* Q2.2.1: constant treatment effect model

* generate variables for polynomials of order 3-6
gen p2 = pov^2
gen p3 = pov^3
gen p4 = pov^4
gen p5 = pov^5
gen p6 = pov^6

* empty matrix to fill with point estimates and standard errors
matrix tbl = J(2,4,.)

* order 3 (run reg, grab ests, predict vals, plot, save for combine later)
reg rel_post t pov p2 p3, vce(hc2)
matrix tbl[1,1] = _b["t"]
matrix tbl[2,1] = _se["t"]
capture drop pred
predict pred
twoway scatter pred pov, title("Order 3")
graph save "s\2_1a.gph", replace

* order 4
reg rel_post t pov p2 p3 p4, vce(hc2)
matrix tbl[1,2] = _b["t"]
matrix tbl[2,2] = _se["t"]
capture drop pred
predict pred
twoway scatter pred pov, title("Order 4")
graph save "s\2_1b.gph", replace

* order 5
reg rel_post t pov p2 p3 p4 p5, vce(hc2)
matrix tbl[1,3] = _b["t"]
matrix tbl[2,3] = _se["t"]
capture drop pred
predict pred
twoway scatter pred pov, title("Order 5")
graph save "s\2_1c.gph", replace

* order 6
reg rel_post t pov p2 p3 p4 p5 p6, vce(hc2)
matrix tbl[1,4] = _b["t"]
matrix tbl[2,4] = _se["t"]
capture drop pred
predict pred
twoway scatter pred pov, title("Order 6")
graph save "s\2_1d.gph", replace

* combine graphs
graph combine "s\2_1a.gph" "s\2_1b.gph" "s\2_1c.gph" "s\2_1d.gph"
graph export "s\2_1.png", replace

* write table for LaTeX
mat2txt, matrix(tbl) saving("s\2_1.txt") format(%9.4f) replace

********************************************************************************
* Q2.2.2: heterogeneous treatment effect model

********************************************************************************
* Q2.2.3: local parametric model

********************************************************************************
* Q2.3.1: MSE-optimal RD estimators

********************************************************************************
* Q2.3.2: Robustness checks

********************************************************************************
* Q2.4: Local Randomization Methods
	
********************************************************************************
log close
********************************************************************************
