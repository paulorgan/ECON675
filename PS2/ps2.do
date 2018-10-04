********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS2
* Last Update: Oct 3, 2018
********************************************************************************
clear all
set more off
capture log close

cd "C:\Users\prorgan\Box\Classes\Econ 675\Problem Sets\PS2"
log using ps2.log, replace
set seed 22

********************************************************************************
*** Question 1: Kernel Density Estimation
********************************************************************************
* Q1.3a

set obs 1000

gen x = .5*rnormal(-1.5,1.5) + .5*rnormal(1,1)

kdensity x, at(x) generate(fx) nograph

forvalues i = 1/1000 {
	kdensity x if _n != `i', at(x) generate(foo fx`i') kernel(epan2) nograph
	drop foo
}

********************************************************************************
