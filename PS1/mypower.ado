********************************************************************************
** Econ-675: Applied Microeconometrics, Fall 2015
** Lecture 1: Causal Inference and Randomized Experiments
** Author: Matias D. Cattaneo
** Last update: 05-Sep-2015
********************************************************************************
** Power Calculation ADO command
********************************************************************************
program drop _all
capture program drop mypower

program mypower, rclass
	set more off	
	syntax varlist(min=2 max=2) [if] [in] [, tau0(real 0) tau(real 1) alpha(real 0.05) ///
											 plot(real 1) rangeL(real -10) rangeU(real 10)]

	tokenize "`varlist'"
	local y = "`1'"
	local t = "`2'"
	
	qui sum `y' if `t'==0
	local N0 = r(N)
	local mu0 = r(mean)
	local sd0 = r(sd)

	qui sum `y' if `t'==1
	local N1 = r(N)
	local mu1 = r(mean)
	local sd1 = r(sd)
	
	local popt = 1/(1+`sd0'/`sd1')

	local se = sqrt(`sd1'^2 / `N1' + `sd0'^2 / `N0' )
	
	local power = 1 - normal((`tau0'-`tau')/`se'+invnormal(1-`alpha'/2)) ///
	                + normal((`tau0'-`tau')/`se'-invnormal(1-`alpha'/2))

	if (`plot' == 1 ) {
		local rangeL = `rangeL'-`tau0'
		local rangeU = `rangeU'-`tau0'
		twoway function y = 1 - normal((`tau0'- x)/`se'+invnormal(1-`alpha'/2)) ///
			                  + normal((`tau0'- x)/`se'-invnormal(1-`alpha'/2)), ///
		range(`rangeL' `rangeU') yline(`alpha', lpattern(dash))
	}
	local powerr = round(`power', .001)
	
	display as text "Optimal number obs is `popt' (H0:tau0=`tau0', alpha=`alpha')."
	display as text "Power for tau = `tau' is `power' (H0:tau0=`tau0', alpha=`alpha')."
	return local power = `power'
end

