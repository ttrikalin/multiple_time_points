// Perform the  conditional analysis of "moving the baseline"
//         This needs the ncymat and ncSmat matrices of the original analyses. 

use data, clear 

drop in 17   // exclude study 17 that has missing values -- to match the
	     // true conditional analyses 


local m = 4   // time points 
local n = 16  // studies -- we exluded study 17
local m1 = 2  // condition on the first m1 timepoints -- in our example
              // m1 = 1 ==> 6 months
              // m1 = 2 ==> 12 months
              // m1 = 3 ==> 18 months

local m2 = `m' - `m1'  // number of effects of interest

// move the baseline
// number alive at the assessment of the m1-th time point 
// are the new total ("move the baseline")
replace N1 = n1`m1'
replace N2 = n2`m1'

forval j=1/`m2' {
	replace n1`j' = n1`=`j'+`m1''
	replace n2`j' = n2`=`j'+`m1''
}

// clean up old variables. 
forval j= 1/`m1' {
	drop n1`=`m'-`j'+1' n2`=`m'-`j'+1' 
}


// effect sizes (log odds ratios) for conditional b's
forval j = 1/`m2' {
	gen ncb`j' = logit(n1`j'/N1) - logit(n2`j'/N2)  
	// continuity correction 
	replace ncb`j' = logit((n1`j'+.5)/(N1+1)) - logit((n2`j'+.5)/(N2+1)) if ncb`j' ==.   

	forval k=`j'/`m2' {
		// variances 
		if (`j'==`k') {
			gen ncV`j'`j' = 1/n1`j' + 1/(N1-n1`j') + 1/n2`j' + 1/(N2-n2`j') 
			
			// continuity correction 
			replace ncV`j'`j' = 1/(n1`j'+.5) + 1/(N1-n1`j'+.5) + ///
					    1/(n2`j'+.5) + 1/(N2-n2`j'+.5) if ncV`j'`j' ==.

		}
		// covariances
		if (`j'<`k') {
			gen ncV`j'`k' = N1/(n1`j' * (N1 - n1`k')) + N2/(n2`j' * (N2-n2`k')) 

			// continuity correction 
			replace ncV`j'`k' = (N1+1)/((n1`j'+0.5) * (N1 - n1`k' + 0.5)) + ///
					    (N2+1)/((n2`j'+0.5) * (N2 - n2`k' + 0.5)) if ncV`j'`k' ==.
		}
	}
}


forval i=1/`n' {
	mat ncymat`i' = J(1, `m2', 0)	
	forval j =1/`m2' {
		qui summ b`j' in `i' , meanonly		 
		mat ncymat`i'[1,`j'] = r(mean) 
	}
	if (matmissing(ncymat`i')) {
		noi di in yellow "Study `i' has missing values:"
		noi mat li ncymat`i'
	}

	mat ncSmat`i' = J(`m2', `m2' , 0)
	forval k =1/`m2' {
		forval l =`k'/`m2' {
			qui summ  V`k'`l' in `i' , meanonly 
			mat ncSmat`i'[`k', `l'] = r(mean)
			mat ncSmat`i'[`l', `k'] = r(mean)
		}
	}	
}


// We must assert that these matrices are positive definite; the fastest way
// would be to try-catch a Cholesky decomposition (See G Strang Linear Algebra)
// but here i will simply check the eigenvalues, as described in the paper. 
// correctmemat is a short program that does this correction quickly.

forval i=1/`n' {
	correctmemat, matname(ncSmat`i') epsilon(0.08)
	if (r(corrected)) mat ncSmat`i' = r(M)
}


// This is the fixed effects model, fit with GLS. The same can be fit using ML
// but I prefer to show this coding also. So bear with me. 

if (0==0) {  
	mat Wsum = J(`m2', `m2', 0)
        mat Wysum = J(1, `m2', 0)

        forval k=1/`n' {
                mat W`k' = syminv(ncSmat`k')
                mat Wysum = Wysum + ncymat`k' * W`k'
                mat Wsum = Wsum + W`k'
        }
        mat condcovbetaFixed = syminv(Wsum)
        mat condbetaFixed = Wysum * condcovbetaFixed
}



// Report the results. 

tempname cond
postfile  `cond' conditional max_q sorter logor se str50 ( timepoint fix_ran uni_multi ) ///
	using ConditionalResults_`=6*`m1''mo,  replace

local k = 0

forval i = 1/`m2' {
	local k = `k' +1

	// the contrast matrix 
	mat a = J(1, `m2', 0)
	mat a[1,`i'] = 1

	// Naive conditional: the effect size and its variance
	mat ES  = a*condbetaFixed'
	mat VAR  = a*(condcovbetaFixed)*a'
	post `cond' (1) (`m2')  (`k') (`=el(ES, 1, 1)') ///
		(`=sqrt(el(VAR,1,1))')  ("`=(`i'+`m1')*6' months, moving baseline to `=6*`m1'' mo") ///
		("fixed") ("multivariate")
	
}

postclose `cond'


use ConditionalResults_`=6*`m1''mo, clear

sort conditional , stable

// significance level
local alpha = 0.05

// calculate the uncorrected CIs
gen l_uncorrected = logor - abs(invnorm(`alpha'/2)) * se 
gen h_uncorrected = logor + abs(invnorm(`alpha'/2)) * se 

// For multivariate use simutaneous CI's 
gen l_corrected = logor - sqrt(invchi2(`m2', 1-`alpha')) * se 
gen h_corrected = logor + sqrt(invchi2(`m2', 1-`alpha')) * se


gen ORCorrStr = string(exp(logor), "%9.2f") + " (" + string(exp(l_corrected), "%9.2f") + ///
	", " + string(exp(h_corrected), "%9.2f") + ")"
gen ESCorrStr = string((logor), "%9.2f") + " (" + string((l_corrected), "%9.2f") + ///
	", " + string((h_corrected), "%9.2f") + ")"

gen ORUncorrStr = string(exp(logor), "%9.2f") + " (" + string(exp(l_uncorrected), "%9.2f") + ///
	", " + string(exp(h_uncorrected), "%9.2f") + ")"
gen ESUncorrStr = string((logor), "%9.2f") + " (" + string((l_uncorrected), "%9.2f") + ///
	", " + string((h_uncorrected), "%9.2f") + ")"

compress

save ConditionalResults_`=6*`m1''mo, replace 

