// Stata code for main analyses described in Trikalinos and Olkin
// Meta-analysis of outcomes at multiple follow up times: a multivariate approach
// 

// here we demonstrate how to perform conditional analyses.
// as discussed in the paper I am first calculating the conditional 
// covariance matrices for the logit of the proportion of survivors in treatments and controls
// For simplicity of exposition we condition on having survived 6 months

// Because study 17 has missing values for 6 and 18 months we cannot use it in 
// the conditional analysis -- short of imputing 
// So here I exclude it and proceed with 16 studies:

use data, clear 
local n = 16
local p = 4


// Step 1. Form the vectors of the transformed proportions and their covariance matrices 
//         in the treatment and control arms 

forval j=1/`p' {
	gen t`j' = logit(n1`j'/N1)
	replace t`j' = logit((n1`j'+0.5)/(N1+1)) if t`j' == .  // continuity correction

	gen c`j' = logit(n2`j'/N2)
	replace c`j' = logit((n2`j'+0.5)/(N2+1)) if c`j' == .  // continuity correction

	forval k=`j'/`p' {
		if (`j' == `k') {   // variances 
			gen tV`j'`k' = 1/n1`j' + 1/(N1-n1`j')
			replace tV`j'`k' = 1/(n1`j'+0.5) + 1/(N1-n1`j'+0.5) if tV`j'`k'==.

			gen cV`j'`k' = 1/n2`j' + 1/(N2-n2`j')
			replace cV`j'`k' = 1/(n2`j'+0.5) + 1/(N2-n2`j'+0.5) if cV`j'`k'==.
		}
		if (`j'<`k') {   // covariances (keep them separate from the variances)
			gen tV`j'`k' = N1/(n1`j'*(N1-n1`k'))
			replace tV`j'`k' = (N1+1)/((n1`j'+0.5)*(N1-n1`k'+0.5)) if tV`j'`k'==.

			gen cV`j'`k' = N2/(n2`j'*(N2-n2`k'))
			replace cV`j'`k' = (N2+1)/((n2`j'+0.5)*(N2-n2`k'+0.5)) if cV`j'`k'==.
		}
	}
}



forval i =1/`n' {
	mat tVec`i' = J(`p', 1, 0)
	mat tCovMat`i' = J(`p', `p', 0)

	mat cVec`i' = J(`p', 1,  0)
	mat cCovMat`i' = J(`p', `p', 0)

	forval j=1/`p' {
		mat tVec`i'[`j',1] = `=t`j'[`i']' 
		mat cVec`i'[`j',1] = `=c`j'[`i']'

		forval k=`j'/`p' {
			if (`j'<=`k') {
				mat tCovMat`i'[`j',`k'] =  `=tV`j'`k'[`i']'
				mat tCovMat`i'[`k',`j'] = tCovMat`i'[`j', `k'] 

				mat cCovMat`i'[`j',`k'] =  `=cV`j'`k'[`i']'
				mat cCovMat`i'[`k',`j'] = cCovMat`i'[`j', `k'] 
			}
		}
	}
}

// Step 2. Assert that the resulting covariance matrices are positive definite
//         This uses the same method we discussed in the main analyses. 
//	   Choice of epsilon should be examined with plots approach
//         but here we only make an exposition 
foreach X in t c {
	forval i=1/`n' {
		correctmemat , matname(`X'CovMat`i') epsilon(1e-3)
		mat `X'CovMat`i' = r(M)
	}
}


// Step 3. Create the conditional means and covariance matrices as described in the paper
//         If necessary, use a permutation matrix to move the q variables on which we condition
//         first, followed by the rest. Her we need not do that, we condition on 6 months 
//         which is the first time point. 

local q = 1 // we condition on the q=1 st time point.   

forval i = 1/`n' {
	foreach X in t c {
		// these are the partitions of the XCovMat
		// we overwrite the matrices all the time 
		mat Sigma_21 = `X'CovMat`i'[`=`q'+1'..`p', 1..`q']
		mat Sigma_11 = `X'CovMat`i'[1..`q', 1..`q']
		mat Sigma_12 = `X'CovMat`i'[1..`q', `=`q'+1'..`p']
		mat Sigma_22 = `X'CovMat`i'[`=`q'+1'..`p', `=`q'+1'..`p']
		
		// the conditional means are P_2.1 = P_2 - Sigma_21 * inv(Sigma_11) * P_1
		mat `X'CondMeans`i' = `X'Vec`i'[`=`q'+1'..`p',1] - ///
				Sigma_21 * syminv(Sigma_11) * `X'Vec`i'[1..`q',1]
				
	
		// the conditional covariance matrix is 
		// Sigma_22.1 = Sigma_22 - Sigma_21 *inv(Sigma_11)*Sigma_12
	
		mat `X'CondCov`i' = Sigma_22 - Sigma_21 * syminv(Sigma_11) * Sigma_12
	}
}

// Step 4. Form the differences vector and the covariance matrices 
forval i=1/`n' {
	// note -- make the differences a row vector for Step 5...
	mat condDiff`i' = tCondMeans`i'' - cCondMeans`i''  
	mat condCov`i' = tCondCov`i' + cCondCov`i' 
}

// Save the conditional vectors in the dataset 
forval j=1/`=`p'-`q'' {
	gen cb`j' = .a
	forval i=1/`n' {
		replace cb`j' = el(condDiff`i', 1, `j') in `i'
	}
}


// Step 5. repeat the multivariate analyses in the same way

// This is the fixed effects model, fit with GLS. The same can be fit using REML
// but I prefer to show this coding also. So bear with me. 
// Because i have no covariates, there is no design matrix here (it's I(4), omitted)

// now we are left with m= p-q time points (we conditioned on 6 months)
local m = `p'-`q'

if (0==0) {  
	mat Wsum = J(`m', `m', 0)
        mat Wysum = J(1, `m', 0)

        forval k=1/`n' {
                mat W`k' = syminv(condCov`k')
                mat Wysum = Wysum + condDiff`k' * W`k'
                mat Wsum = Wsum + W`k'
        }
        mat condcovbetaFixed = syminv(Wsum)
        mat condbetaFixed = Wysum * condcovbetaFixed
}



// Report the results. 
// remember m= p-q


tempname cond
postfile  `cond' naive max_q sorter logor se str50 ( timepoint fix_ran uni_multi ) ///
	using ConditionalResults_`=6*`q''mo,  replace

local k = 0

forval i = `=`q'+1'/`p' {
	local k = `k' +1

	// the contrast matrix 
	mat a = J(1, `m', 0)
	mat a[1, `=`i'-`q''] = 1

	// Conditional analyses: the effect size and its variance
	mat ES  = a*condbetaFixed'
	mat VAR  = a*(condcovbetaFixed)*a'
	post `cond' (0) (`q') (`k') (`=el(ES, 1, 1)') ///
		(`=sqrt(el(VAR,1,1))')  ("`=`i'*6' months | `=6*(`q'+1)' mo") ///
		("fixed") ("multivariate")
 	
	// Naive conditional: the effect size and its variance
	mat ES  = a*naivecondbetaFixed'
	mat VAR  = a*(naivecondcovbetaFixed)*a'
	post `cond' (1) (`q')  (`k') (`=el(ES, 1, 1)') ///
		(`=sqrt(el(VAR,1,1))')  ("`=`i'*6' months, moving baseline to `=6*(`q'+1)' mo") ///
		("fixed") ("multivariate")
	
}

postclose `cond'


use ConditionalResults_`=6*`q''mo, clear

sort naive , stable

// number of timepoints (outcomes)
// remember m=p-q 

// significance level
local alpha = 0.05

// calculate the uncorrected CIs
gen l_uncorrected = logor - abs(invnorm(`alpha'/2)) * se 
gen h_uncorrected = logor + abs(invnorm(`alpha'/2)) * se 

// For multivariate use simutaneous CI's 
gen l_corrected = logor - sqrt(invchi2(`m', 1-`alpha')) * se 
gen h_corrected = logor + sqrt(invchi2(`m', 1-`alpha')) * se


gen ORCorrStr = string(exp(logor), "%9.2f") + " (" + string(exp(l_corrected), "%9.2f") + ///
	", " + string(exp(h_corrected), "%9.2f") + ")"
gen ESCorrStr = string((logor), "%9.2f") + " (" + string((l_corrected), "%9.2f") + ///
	", " + string((h_corrected), "%9.2f") + ")"

gen ORUncorrStr = string(exp(logor), "%9.2f") + " (" + string(exp(l_uncorrected), "%9.2f") + ///
	", " + string(exp(h_uncorrected), "%9.2f") + ")"
gen ESUncorrStr = string((logor), "%9.2f") + " (" + string((l_uncorrected), "%9.2f") + ///
	", " + string((h_uncorrected), "%9.2f") + ")"

compress

save ConditionalResults_`=6*`q''mo, replace 




