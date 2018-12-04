********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS6
* Last Update: Dec 4, 2018
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
    graph_options(title("Evenly-spaced binning, IMSE-optimal", size(medium)) ///
					legend(size(vsmall)))
graph save "s\1_1a.gph", replace

rdplot rel_pre pov, c(0) binselect(esmv) ///
    graph_options(title("Evenly-spaced binning, variability-mimicking", size(medium)) ///
					legend(size(vsmall)))
graph save "s\1_1b.gph", replace

rdplot rel_pre pov, c(0) binselect(qs) ///
    graph_options(title("Quantile-spaced binning, IMSE-optimal", size(medium)) ///
					legend(size(vsmall)))
graph save "s\1_1c.gph", replace

rdplot rel_pre pov, c(0) binselect(qsmv) ///
    graph_options(title("Quantile-spaced binning, variability mimicking", size(medium)) ///
					legend(size(vsmall)))
graph save "s\1_1d.gph", replace

graph combine "s\1_1a.gph" "s\1_1b.gph" "s\1_1c.gph" "s\1_1d.gph"
graph export "s\1_1s.png", replace

********************************************************************************
* Q2.1.2: Falsification Tests
* falsification of RD design using three methods

gen t = pov > 0

* (1) histogram plots
twoway (hist pov if t, freq bcolor(red)) ///
	   (hist pov if !t, freq bcolor(black)), ///
	   legend(label(1 "Treated") label(2 "Untreated") order(2 1))
graph export "s\1_2s.png", replace
	   
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
graph export "s\2_1s.png", replace

* write table for LaTeX
mat2txt, matrix(tbl) saving("s\2_1.txt") format(%9.4f) replace

********************************************************************************
* Q2.2.2: heterogeneous treatment effect model

* define new, interacted terms for polynomials
gen p1t = pov*t
gen p1u = pov*(1-t)
gen p2t = p2*t
gen p2u = p2*(1-t)
gen p3t = p3*t
gen p3u = p3*(1-t) 
gen p4t = p4*t
gen p4u = p4*(1-t)
gen p5t = p5*t
gen p5u = p5*(1-t)
gen p6t = p6*t
gen p6u = p6*(1-t)

* empty matrix to fill with point estimates and standard errors
matrix tbl = J(2,4,.)

local base = "p1t p1u p2t p2u p3t p3u"

* order 3 (run reg, grab ests, predict vals, plot, save for combine later)
reg rel_post t `base', vce(hc2)
matrix tbl[1,1] = _b["t"]
matrix tbl[2,1] = _se["t"]
capture drop pred
predict pred
twoway scatter pred pov, title("Order 3")
graph save "s\2_2a.gph", replace

* order 4
reg rel_post t `base' p4t p4u, vce(hc2)
matrix tbl[1,2] = _b["t"]
matrix tbl[2,2] = _se["t"]
capture drop pred
predict pred
twoway scatter pred pov, title("Order 4")
graph save "s\2_2b.gph", replace

* order 5
reg rel_post t `base' p4t p4u p5t p5u, vce(hc2)
matrix tbl[1,3] = _b["t"]
matrix tbl[2,3] = _se["t"]
capture drop pred
predict pred
twoway scatter pred pov, title("Order 5")
graph save "s\2_2c.gph", replace

* order 6
reg rel_post t `base' p4t p4u p5t p5u p6t p6u, vce(hc2)
matrix tbl[1,4] = _b["t"]
matrix tbl[2,4] = _se["t"]
capture drop pred
predict pred
twoway scatter pred pov, title("Order 6")
graph save "s\2_2d.gph", replace

* combine graphs
graph combine "s\2_2a.gph" "s\2_2b.gph" "s\2_2c.gph" "s\2_2d.gph"
graph export "s\2_2s.png", replace

* write table for LaTeX
mat2txt, matrix(tbl) saving("s\2_2.txt") format(%9.4f) replace

********************************************************************************
* Q2.2.3: local parametric model (p = 0,1,2 and h = 1,5,9,18)

* matrix to fill with estimates
matrix tbl = J(12,3,.)
matrix tbl[1,1] = 1
matrix tbl[4,1] = 5
matrix tbl[7,1] = 9
matrix tbl[10,1] = 18

* for h = 1
reg rel_post t if abs(pov)<=1, vce(hc2)
matrix tbl[2,1] = _b["t"]
matrix tbl[3,1] = _se["t"]
capture drop pred
predict pred if abs(pov)<=1
twoway scatter pred pov if abs(pov)<=1, title("Order 0, h = 1")
graph save "s\2_3h1a.gph", replace

reg rel_post t p1t p1u if abs(pov)<=1, vce(hc2)
matrix tbl[2,2] = _b["t"]
matrix tbl[3,2] = _se["t"]
capture drop pred
predict pred if abs(pov)<=1
twoway scatter pred pov if abs(pov)<=1, title("Order 1, h = 1")
graph save "s\2_3h1b.gph", replace

reg rel_post t p1t p1u p2t p2u if abs(pov)<=1, vce(hc2)
matrix tbl[2,3] = _b["t"]
matrix tbl[3,3] = _se["t"]
capture drop pred
predict pred if abs(pov)<=1
twoway scatter pred pov if abs(pov)<=1, title("Order 2, h = 1")
graph save "s\2_3h1c.gph", replace

* for h = 5
reg rel_post t if abs(pov)<=5, vce(hc2)
matrix tbl[5,1] = _b["t"]
matrix tbl[6,1] = _se["t"]
capture drop pred
predict pred if abs(pov)<=5
twoway scatter pred pov if abs(pov)<=5, title("Order 0, h = 5")
graph save "s\2_3h5a.gph", replace

reg rel_post t p1t p1u if abs(pov)<=5, vce(hc2)
matrix tbl[5,2] = _b["t"]
matrix tbl[6,2] = _se["t"]
capture drop pred
predict pred if abs(pov)<=5
twoway scatter pred pov if abs(pov)<=5, title("Order 1, h = 5")
graph save "s\2_3h5b.gph", replace

reg rel_post t p1t p1u p2t p2u if abs(pov)<=5, vce(hc2)
matrix tbl[5,3] = _b["t"]
matrix tbl[6,3] = _se["t"]
capture drop pred
predict pred if abs(pov)<=5
twoway scatter pred pov if abs(pov)<=5, title("Order 2, h = 5")
graph save "s\2_3h5c.gph", replace

* for h = 9
reg rel_post t if abs(pov)<=9, vce(hc2)
matrix tbl[8,1] = _b["t"]
matrix tbl[9,1] = _se["t"]
capture drop pred
predict pred if abs(pov)<=9
twoway scatter pred pov if abs(pov)<=9, title("Order 0, h = 9")
graph save "s\2_3h9a.gph", replace

reg rel_post t p1t p1u if abs(pov) <=9, vce(hc2)
matrix tbl[8,2] = _b["t"]
matrix tbl[9,2] = _se["t"]
capture drop pred
predict pred if abs(pov)<=9
twoway scatter pred pov if abs(pov)<=9, title("Order 1, h = 9")
graph save "s\2_3h9b.gph", replace

reg rel_post t p1t p1u p2t p2u if abs(pov)<=9, vce(hc2)
matrix tbl[8,3] = _b["t"]
matrix tbl[9,3] = _se["t"]
capture drop pred
predict pred if abs(pov)<=9
twoway scatter pred pov if abs(pov)<=9, title("Order 2, h = 9")
graph save "s\2_3h9c.gph", replace

* for h = 18
reg rel_post t if abs(pov)<=18, vce(hc2)
matrix tbl[11,1] = _b["t"]
matrix tbl[12,1] = _se["t"]
capture drop pred
predict pred if abs(pov)<=18
twoway scatter pred pov if abs(pov)<=18, title("Order 0, h = 18")
graph save "s\2_3h18a.gph", replace

reg rel_post t p1t p1u if abs(pov)<=18, vce(hc2)
matrix tbl[11,2] = _b["t"]
matrix tbl[12,2] = _se["t"]
capture drop pred
predict pred if abs(pov)<=18
twoway scatter pred pov if abs(pov)<=18, title("Order 1, h = 18")
graph save "s\2_3h18b.gph", replace

reg rel_post t p1t p1u p2t p2u if abs(pov)<=18, vce(hc2)
matrix tbl[11,3] = _b["t"]
matrix tbl[12,3] = _se["t"]
capture drop pred
predict pred if abs(pov)<=18
twoway scatter pred pov if abs(pov)<=18, title("Order 2, h = 18")
graph save "s\2_3h18c.gph", replace

* combine all graphs
graph combine "s\2_3h1a.gph" "s\2_3h1b.gph" "s\2_3h1c.gph" ///
	"s\2_3h5a.gph" "s\2_3h5b.gph" "s\2_3h5c.gph" ///
	"s\2_3h9a.gph" "s\2_3h9b.gph" "s\2_3h9c.gph" ///
	"s\2_3h18a.gph" "s\2_3h18b.gph" "s\2_3h18c.gph", c(3)
graph export "s\2_3s.png", replace

* write table to LaTeX
mat2txt, matrix(tbl) saving("s\2_3.txt") format(%9.4f) replace

********************************************************************************
* Q2.3.1: MSE-optimal RD estimators

* just manually typing these into LaTeX
* 'all' gives us the three methods, p and q tell the polynomial orders
rdrobust rel_post pov, p(0) q(1) c(0) all
rdrobust rel_post pov, p(1) q(2) c(0) all
rdrobust rel_post pov, p(2) q(3) c(0) all

********************************************************************************
* Q2.3.2: Robustness checks

** a) Placebo outcome tests
* check preintervention related mortality
rdrobust rel_pre pov, p(1) q(2) c(0) all
* check postintervention unrelated mortality
rdrobust inj_post pov, p(1) q(2) c(0) all

** b) Bandwidth and Kernel sensitivity
* empty matrix to fill with results
matrix tbl = J(3, 10, .)

* loop over bandwidths, calc for three different kernels
forvalues h = 1(1)10 {
    quiet rdrobust rel_post pov, h(`h') kernel(tri) p(1) q(2) c(0) all
	matrix tbl[1, `h'] = round(e(tau_bc), .01)
	quiet rdrobust rel_post pov, h(`h') kernel(uni) p(1) q(2) c(0) all
	matrix tbl[2, `h'] = round(e(tau_bc), .01)
	quiet rdrobust rel_post pov, h(`h') kernel(epa) p(1) q(2) c(0) all
	matrix tbl[3, `h'] = round(e(tau_bc), .01)
}

* write to LaTeX
mat2txt, matrix(tbl) saving("s\3_2b.txt") format(%9.4f) replace

** c) "donut hole" approach
* first sort by closeness for povrate to 0
gen pov_abs = abs(pov)
sort pov_abs
drop pov_abs

* empty matrix to fill with results
matrix tbl = J(1, 10, .)

* loop over 10 l's
forvalues l = 1(1)10 {
    quiet rdrobust rel_post pov if _n>`l', p(1) q(2) c(0) all
	matrix tbl[1, `l'] = round(e(tau_bc), .01)
}

* write to LaTeX
mat2txt, matrix(tbl) saving("s\3_2c.txt") format(%9.4f) replace

** d) Placebo cutoff approach
* empty matrix to fill with results
matrix tbl = J(2, 10, .)

* loop over various cutoffs, save point estimates and p-values
forvalues c = -10(2)10 {
    if `c' != 0 {
	    quiet rdrobust rel_post pov, c(`c') p(1) q(2) all
		* to know which matrix position to store in
		if `c' < 0 { 
		    local i = `c'/2 + 6 
		}
		if `c' > 0 { 
		    local i = `c'/2 + 5 
		}
		matrix tbl[1, `i'] = round(e(tau_bc), .01)
		matrix tbl[2, `i'] = round(e(pv_rb), .01)
	}
}

* write to LaTeX
mat2txt, matrix(tbl) saving("s\3_2d.txt") format(%9.4f) replace

********************************************************************************
* Q2.4.1: Local Randomization Methods - Window Selection

* select windows using the four included methods
rdwinselect pov inj_post rel_pre, ///
    cutoff(0) wmin(0.1) wstep(0.2) nwindows(20) statistic(diffmeans)
* recommended window is 1.7
	
rdwinselect pov inj_post rel_pre, ///
    cutoff(0) wmin(0.1) wstep(0.2) nwindows(20) statistic(ksmirnov)
* recommended window is 1.9

rdwinselect pov inj_post rel_pre, ///
    cutoff(0) wmin(0.1) wstep(0.2) nwindows(20) statistic(ranksum)
* recommended window is 1.7

rdwinselect pov inj_post rel_pre, ///
    cutoff(0) wmin(0.1) wstep(0.2) nwindows(20) statistic(hotelling)
* recommended window is 1.9

* plot based on the selected binwidth of 1.7 (using diffmeans)
rdplot rel_post pov if abs(pov)<=1.7, c(0) p(0) binselect(es)
graph export "s\4_1s.png", replace

********************************************************************************
* Q2.4.2: Local Randomization Methods - Basic Analysis

* do the randomization analysis using window from above
rdrandinf rel_post pov, wl(-1.7) wr(1.7) seed(22)

********************************************************************************
* Q2.4.2: Local Randomization Methods - Sensitivity Analysis

* empty matrix to fill with results
matrix tbl = J(4, 10, .)

* loop over window values, store results
forvalues i = 1(1)10 {
	local win = 0.8+(`i'-1)*0.2
	matrix tbl[1,`i'] = `win'
	
    reg rel_post t if abs(pov)<=`win', vce(hc2)
	matrix temp = r(table)
    
    matrix tbl[2, `i'] = temp["b","t"]
    matrix tbl[3, `i'] = temp["se","t"]
    matrix tbl[4, `i'] = temp["pvalue","t"]
}

* write to LaTeX
mat2txt, matrix(tbl) saving("s\4_3.txt") format(%9.4f) replace
	
********************************************************************************
log close
********************************************************************************
