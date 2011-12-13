// Here we compare univariate and multivariate fixed effects meta-analysis

use data, clear
stack b1 V11 b2 V22 b3 V33 b4 V44 , into(b V) clear
gen se = sqrt(V)

tempname uni 

postfile  `uni' sorter logor se Q df str20 ( timepoint fix_ran uni_multi) using Results, replace

forval i = 1/4 {
	qui metan b se if _stack ==`i' , fixed nograph  
	post `uni' (`i') (r(ES)) (r(seES)) (r(het)) (r(df)) ///
		("`=6*`i'' months") ("fixed") ("univariate")
	//qui metan b es if _stack ==`i' , random nograph  
	//post `uni' (`i') (r(ES)) (r(seES)) (r(het)) (r(df)) ///
	//	("`=6*`i'' months") ("random") ("univariate")
    qui metareg b if _stack ==`i' , wsse(se) z 
    mat bb= e(b)
    mat vv= e(V)
    post `uni' (`i') (`=el(bb,1,1)') (`=sqrt(el(vv,1,1))') (e(Q)) (e(df_Q)) ///
		("`=6*`i'' months") ("random") ("univariate")
}

postclose `uni'

