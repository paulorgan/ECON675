********************************************************************************
* Author: Paul R. Organ
* Purpose: ECON 675, PS2
* Last Update: Oct 11, 2018
*
* NOTE: For PS2 I mostly focused on implementation in R
* This code relies heavily on my classmates' work
********************************************************************************
clear all
set more off
set matsize 10000
capture log close

cd "C:\Users\prorgan\Box\Classes\Econ 675\Problem Sets\PS2"
log using ps2.log, replace

********************************************************************************
*** Question 1: Kernel Density Estimation
********************************************************************************

* (h values from R)
global n = 1000 
global hvalues 0.4099711 0.4919653 0.5739595 0.6559538 0.7379480 0.8199422 0.9019364 0.9839307 1.0659249 1.1479191 1.2299133
global nh = 11
global M = 1000 //number of iterations
mat hvalues = (0.4099711, 0.4919653, 0.5739595, 0.6559538, 0.7379480, 0.8199422, 0.9019364, 0.9839307, 1.0659249, 1.1479191, 1.2299133)

* Caluclate MSE
	mata: 
	
	//**********FUNCTIONS***********
	// function for calculating kernel 
	real scalar function kern(real scalar u){
		return(.75*(1-u^2)*(abs(u)<=1))
		}
	
	// function for calculating true density
	real scalar function f_true(real scalar u){
		return(.5*normalden(u,-1.5,sqrt(1.5)) + .5*normalden(u,1,1))
		}	
	
	// function for calculating MSE (LI & LO)
	real vector function mse(real vector xdata, real scalar hvalue){
		//Construct two matrices of xdata
		M1 = J($n,$n,.) // n x n matrix with one column for each observation
		M2 = J($n,$n,.) // n x n matrix with one row for each observation
		for (i=1; i<= $n; i++) {
			v = J($n,1,xdata[i])
			M1[,i] = v
			M2[i,] = v'
			}
		
		M3 = (M1-M2)/hvalue //object to be evaluated by kernel
		M4 = J($n,$n,.)
		M5 = J($n,$n,.)
		fx = J($n,1,.)
		
		for (i=1; i<=$n; i++){
			for (j=1; j<=$n; j++){
				M4[i,j] = kern(M3[i,j])
			}
			M5[i,] = M4[i,]
			M5[i,i]=0
			
			fx[i,1] = f_true(xdata[i])
		}
		
		fhat_LI = rowsum(M4)/($n*hvalue)
		fhat_LO = rowsum(M5)/(($n-1)*hvalue)
		
		sqe_LI = (fhat_LI-fx):^2
		sqe_LO = (fhat_LO-fx):^2
		
		mse_LI = mean(sqe_LI)
		mse_LO = mean(sqe_LO)
		
		return((mse_LI,mse_LO))
		}
	
	// function for importing/exporting to mata for mse calculation
	void iteration(real scalar m){
		x= st_data(.,.)
		hvalues = st_matrix("hvalues")
		
		mse = J(11,2,.)
		for (h=1; h<=11; h++){
			mse[h,] = mse(x,hvalues[1,h])
		}	
		st_matrix("msetemp",mse)
		}
	end

* DGP
	*Normal density specs
	global mu1 = -1.5
	global mu2 = 1
	global sd1 = sqrt(1.5)
	global sd2 = 1
	
	*Empty matrix to be filled
	mat msesum = J(11,2,0)
	
	*Loop through iterations
	timer on 1
	forval m = 1/$M{
	disp `m'
	set obs $n
	
	*equally weight two normal distributions
	gen comps = uniform() >= .5

	*generate sample 
	gen x = comps*rnormal($mu1,$sd1) + (1-comps)*rnormal($mu2,$sd2)
	drop comps
	
	*call mata function to calculate mse
	mata iteration(`m')
	mat msesum = msesum + msetemp
	drop x
	}
timer off 1
timer list

mat imse = msesum/1000
svmat imse
rename imse1 imse_li
rename imse2 imse_lo

egen h = fill(0.4099711 0.4919653 0.5739595 0.6559538 0.7379480 0.8199422 0.9019364 0.9839307 1.0659249 1.1479191 1.2299133)

twoway(line imse_li h)(line imse_lo h)
graph export q1_3b_S.png, replace

********************************************************************************
*** Question 2: Linear Smoothers, Cross-Validation, and Series
********************************************************************************
* clean up
clear all

********************************************************************************
* Q2.5b
set obs 1000

* mata function to calculate CV statistic
* CV(list, i): list=variable list, i = max polynomial
mata
void CV(list, i) {
st_view(y=., ., "y")
st_view(X=., ., tokens(list))
XpX  = cross(X, X)
XpXinv  = invsym(XpX)
b  = XpXinv*cross(X, y)
w = diagonal(X*XpXinv*X')
muhat = X*b
num = (y - muhat):*(y - muhat)
den= (J(1000,1,1) - w):*(J(1000,1,1) - w)
div = num:/den
CV = mean(div)
CV
st_numscalar("mCV"+strofreal(i), CV)
}
end

* Program which runs the monte-carlo experiment
program CVsim, rclass
drop _all
set obs 1000
forvalues i = 0/20 { 
gen CV`i' = 0
}
gen x = runiform(-1,1)
gen e = x^2*(rchi2(5)-5)
gen y = exp(-0.1*(4*x-1)^2)*sin(5*x)+e
forvalues i = 0/20 {
gen x`i' = x^`i'
}
forvalues i = 0/20 {
	global xlist = "x0-x`i'"
	di "$xlist"
	mata CV("$xlist", `i')
	replace CV`i' = mCV`i'
	}
end 

* Run the experiment
set seed 22
simulate CV0=CV0 CV1=CV1 CV2=CV2 CV3=CV3 CV4=CV4 CV5=CV5 CV6=CV6 CV7=CV7 CV8=CV8 /// 
		 CV9=CV9 CV10=CV10 CV11=CV11 CV12=CV12 CV13=CV13 CV14=CV14 CV15=CV15 ///
		 CV16=CV16 CV17=CV17 CV18=CV18 CV19=CV19 CV20=CV20, reps(100) nodots: CVsim
collapse *
gen i = 1

reshape long CV, i(i) j(k)
sort CV
local min = k[1]
twoway scatter CV k, ytitle("Mean CV") xtitle("K") xlabel(0(2)20) xmtick(0(1)20) xline(`min')

graph export q2_5b_S.png, replace

********************************************************************************
* Q2.5c
* Program which runs the monte-carlo experiment. At the beginning we generate new
* data, then we call the mata command
program muhatsim, rclass
drop _all
set obs 1000
gen x = runiform(-1,1)
gen e = x^2*(rchi2(5)-5)
gen y = exp(-0.1*(4*x-1)^2)*sin(5*x)+e
forvalues p = 0/7 { 
gen x`p' = x^`p'
}
reg y x0-x7, nocons
clear
set obs 11
gen n = _n
gen foo = 1
gen x = -1+(_n-1)/5
forvalues p = 0/7 { 
gen x`p' = x^`p'
}
predict muhat
predict se, stdp
generate lb = muhat - invnormal(0.975)*se
generate ub = muhat + invnormal(0.975)*se
// gen muhat = _b[x0]*xgp0 + _b[x1]*xgp1 + _b[x2]*xgp2+ _b[x3]*xgp3 + _b[x4]*xgp4	///
//			+ _b[x5]*xgp5 + _b[x6]*xgp6 + _b[x7]*xgp7
keep n muhat foo lb ub 
reshape wide muhat lb ub, i(foo) j(n)
end
set seed 22
simulate muhat1=muhat1 muhat2=muhat2 muhat3=muhat3 muhat4=muhat4 muhat5=muhat5 ///
	muhat6=muhat6 muhat7=muhat7 muhat8=muhat8 muhat9=muhat9 muhat10=muhat10 muhat11=muhat11 ///
	ub1=ub1 ub2=ub2 ub3=ub3 ub4=ub4 ub5=ub5 ub6=ub6 ub7=ub7 ub8=ub8 ub9=ub9 ub10=ub10 ub11=ub11 ///
	lb1=lb1 lb2=lb2 lb3=lb3 lb4=lb4 lb5=lb5 lb6=lb6 lb7=lb7 lb8=lb8 lb9=lb9 lb10=lb10 lb11=lb11, reps(1000) nodots: muhatsim
gen i = _n
reshape long muhat ub lb, i(i) j(grid)
collapse muhat ub lb, by(grid)
gen x = -1+ (grid-1)/5
twoway (function y = exp(-0.1*(4*x-1)^2)*sin(5*x), range(-1 1) lcolor(red)) ///
		(line muhat x, lcolor(gs6)) (line lb x, lcolor(gs6) lpattern(dash)) (line ub x, lcolor(gs6) lpattern(dash)), ///
		legend(order(1 "True f" 2 "Estimated f" 3 "Confidence Interval") rows(1)) ytitle(Y) xtitle(X)
graph export q2_5c_S.png, replace
********************************************************************************
