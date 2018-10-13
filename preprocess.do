/*
	First, create a set with each decision as a row in the data, all
	information about individuals stripped out.
*/
use data/DataErik, clear
drop session round spectator insurance risk gender 
reshape long dist, i(pid) j(situation)
merge m:1 situation using data/situations
drop _merge
label var dist "Equalize (indicator)"
gen int inc1 = (dist==1)*X/2 + (dist==0)*e1
gen int inc2 = (dist==1)*X/2 + (dist==0)*e2
drop dist
label var inc1 "Income of P1"
label var inc2 "Income of P2"
order pid situation outcome1 outcome2 e1 e2 inc1 inc2 X Fse Fl Fle Fscc
compress
sort pid situation
saveold data/DataErik_long, replace version(12)

/*
	Second, create a dataset with info about participants, which session
	and such.
*/
use data/DataErik, clear
keep pid session round spectator insurance risk gender
encode session, g(esession)
drop session 
rename esession Session
gen byte sex = 1 + (gender=="Female")
label define sex 1 "Male" 2 "Female"
label values sex sex
drop gender
order pid Session round sex spectator insurance risk
compress
sort pid
saveold data/DataErik_players, replace version(12)


