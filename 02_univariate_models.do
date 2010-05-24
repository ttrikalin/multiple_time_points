// Stata code for main analyses described in Trikalinos and Olkin
// Meta-analysis of outcomes at multiple follow up times: a multivariate approach
// 

// Here we compare univariate iand multivariate fixed effects meta-analysis
// You will need the metan routine in Stata

use data, clear
stack b1 V11 b2 V22 b3 V33 b4 V44 , into(b V) clear
gen es = sqrt(V)

tempname uni 

postfile  `uni' sorter logor se Q df str20 ( timepoint fix_ran uni_multi) using Results, replace

forval i = 1/4 {
	qui metan b es if _stack ==`i' , fixed nograph  
	post `uni' (`i') (r(ES)) (r(seES)) (r(het)) (r(df)) ///
		("`=6*`i'' months") ("fixed") ("univariate")
	qui metan b es if _stack ==`i' , random nograph  
	post `uni' (`i') (r(ES)) (r(seES)) (r(het)) (r(df)) ///
		("`=6*`i'' months") ("random") ("univariate")
}

postclose `uni'

