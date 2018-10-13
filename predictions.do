set more off

/*
	Parameters for simulation:
*/ 
set seed 123456789
global NS=100000
global l_se = 0.343
global l_l = 0.314
global l_cc = 0.343
global mu=-0.912
global sigma =2.140


/*
	Simulating individual parameters:
*/
use data/DataErik_players, clear
keep if round==1
expand ${NS}
tempvar subpid
bys pid: gen `subpid' = _n
gen long spid = pid*1000000 + `subpid'
gen double u = runiform()
gen byte p_se = (u<${l_se})
gen byte p_cc = (u>= (1- ${l_cc}))
gen byte p_l  = (p_se==0 & p_cc==0)
gen double gamma = exp( ${mu} + ${sigma}*rnormal()) 
keep spid p_se p_l p_cc gamma
tempfile params
save `params'

/*
	Simulating choices:
*/
use data/DataErik_long, clear

expand ${NS}
bys situation pid: gen s1 = _n
gen long spid = pid*1000000 + s1
merge m:1 spid using `params'
keep if _merge==3
drop _merge
gen e_equal = -log(-log(runiform()))
gen e_leave = -log(-log(runiform()))
gen    ud_equal = (-(gamma /X)*(X/2 - Fse)^2 + e_equal) - (-(gamma/X)*(e1-Fse)^2 + e_leave) if p_se==1
replace ud_equal= (-(gamma /X)*(X/2 - Fl )^2 + e_equal) - (-(gamma/X)*(e1-Fl )^2 + e_leave) if p_l ==1
replace ud_equal= (-(gamma /X)*(X/2 - Fcc)^2 + e_equal) - (-(gamma/X)*(e1-Fcc)^2 + e_leave) if p_cc==1
gen equal = (ud_equal>=0)

gen int y1 = X/2 if equal==1
replace y1 = e1 if equal==0
gen cSE = (y1==Fse)
gen cCC = (y1==Fcc)
gen cL  = (y1==Fl)
gen cLE = (y1==Fle)
replace cLE = 1 if inlist(situation,2,5) 
collapse (sum) cSE cLE cCC cL , by(spid) 
forvalues t=8/11 {
	gen cSE`t' = cSE>=`t'
	gen cCC`t' = cCC>=`t'
	gen cLE`t' = cLE>=`t'
	gen cL`t'  = cL >=`t'
}

format %4.3f *11 *10 *9

forvalues t=11(-1)9 {
	di "t=`t'"
	summ cSE`t' cL`t' cLE`t' cCC`t', format
}


