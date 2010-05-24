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
local K = 16  // studies 
local m = 4   // follow up points 


// Step 1. Form the vectors of the transformed proportions and their covariance matrices 
//         in the treatment and control arms 

forval j1=1/`m' {
	gen t`j1' = logit(n1`j1'/N1)
	replace t`j1' = logit((n1`j1'+0.5)/(N1+1)) if t`j1' == .  // continuity correction

	gen c`j1' = logit(n2`j1'/N2)
	replace c`j1' = logit((n2`j1'+0.5)/(N2+1)) if c`j1' == .  // continuity correction

	forval j2=`j1'/`m' {
		if (`j1' == `j2') {   // variances 
			gen tV`j1'`j2' = 1/n1`j1' + 1/(N1-n1`j1')
			replace tV`j1'`j2' = 1/(n1`j1'+0.5) + 1/(N1-n1`j1'+0.5) if tV`j1'`j2'==.

			gen cV`j1'`j2' = 1/n2`j1' + 1/(N2-n2`j1')
			replace cV`j1'`j2' = 1/(n2`j1'+0.5) + 1/(N2-n2`j1'+0.5) if cV`j1'`j2'==.
		}
		if (`j1'<`j2') {   // covariances (keep them separate from the variances)
			gen tV`j1'`j2' = N1/(n1`j1'*(N1-n1`j2'))
			replace tV`j1'`j2' = (N1+1)/((n1`j1'+0.5)*(N1-n1`j2'+0.5)) if tV`j1'`j2'==.

			gen cV`j1'`j2' = N2/(n2`j1'*(N2-n2`j2'))
			replace cV`j1'`j2' = (N2+1)/((n2`j1'+0.5)*(N2-n2`j2'+0.5)) if cV`j1'`j2'==.
		}
	}
}



forval i =1/`K' {
	mat tVec`i' = J(`m', 1, 0)
	mat tCovMat`i' = J(`m', `m', 0)

	mat cVec`i' = J(`m', 1,  0)
	mat cCovMat`i' = J(`m', `m', 0)

	forval j1=1/`m' {
		mat tVec`i'[`j1',1] = `=t`j1'[`i']' 
		mat cVec`i'[`j1',1] = `=c`j1'[`i']'

		forval j2=`j1'/`m' {
			if (`j1'<=`j2') {
				mat tCovMat`i'[`j1',`j2'] =  `=tV`j1'`j2'[`i']'
				mat tCovMat`i'[`j2',`j1'] = tCovMat`i'[`j1', `j2'] 

				mat cCovMat`i'[`j1',`j2'] =  `=cV`j1'`j2'[`i']'
				mat cCovMat`i'[`j2',`j1'] = cCovMat`i'[`j1', `j2'] 
			}
		}
	}
}

// Step 2. Assert that the resulting covariance matrices are positive definite
//         This uses the same method we discussed in the main analyses. 
//	   Choice of epsilon should be examined with plots approach
//         but here we only make an exposition 
foreach X in t c {
	forval i=1/`K' {
		correctmemat , matname(`X'CovMat`i') epsilon(1e-3)
		mat `X'CovMat`i' = r(M)
	}
}


// Step 3. Create the conditional means and covariance matrices as described in the paper
//         If necessary, use a permutation matrix to move the q variables on which we condition
//         first, followed by the rest. Her we need not do that, we condition on 6 months 
//         which is the first time point. 

local m1 = 1 // we condition on the m1=1 st time point.   

forval i = 1/`K' {
	foreach X in t c {
		// these are the partitions of the XCovMat
		// we overwrite the matrices all the time 
		mat Sigma_21 = `X'CovMat`i'[`=`m1'+1'..`m', 1..`m1']
		mat Sigma_11 = `X'CovMat`i'[1..`m1', 1..`m1']
		mat Sigma_12 = `X'CovMat`i'[1..`m1', `=`m1'+1'..`m']
		mat Sigma_22 = `X'CovMat`i'[`=`m1'+1'..`m', `=`m1'+1'..`m']
		
		// the conditional means are P_2.1 = P_2 - Sigma_21 * inv(Sigma_11) * P_1
		mat `X'CondMeans`i' = `X'Vec`i'[`=`m1'+1'..`m',1] - ///
				Sigma_21 * syminv(Sigma_11) * `X'Vec`i'[1..`m1',1]
				
	
		// the conditional covariance matrix is 
		// Sigma_22.1 = Sigma_22 - Sigma_21 *inv(Sigma_11)*Sigma_12
	
		mat `X'CondCov`i' = Sigma_22 - Sigma_21 * syminv(Sigma_11) * Sigma_12
	}
}

// Step 4. Form the differences vector and the covariance matrices 
forval i=1/`K' {
	// note -- make the differences a row vector for Step 5...
	mat condDiff`i' = tCondMeans`i'' - cCondMeans`i''  
	mat condCov`i' = tCondCov`i' + cCondCov`i' 
}

// Save the conditional vectors in the dataset 
forval j1=1/`=`m'-`m1'' {
	gen cb`j1' = .a
	forval j2=1/`K' {
		replace cb`j1' = el(condDiff`j2', 1, `j1') in `j2'
	}
}


// Step 5. repeat the multivariate analyses in the same way

// This is the fixed effects model, fit with GLS. The same can be fit using REML
// but I prefer to show this coding also. So bear with me. 
// Because i have no covariates, there is no design matrix here (it's I(4), omitted)

// now we are left with m2= m-m1 time points (we conditioned on 6 months)
local m2 = `m'-`m1'

if (0==0) {  
	mat Wsum = J(`m2', `m2', 0)
        mat Wysum = J(1, `m2', 0)

        forval j=1/`K' {
                mat W`j' = syminv(condCov`j')
                mat Wysum = Wysum + condDiff`j' * W`j'
                mat Wsum = Wsum + W`j'
        }
        mat condcovbetaFixed = syminv(Wsum)
        mat condbetaFixed = Wysum * condcovbetaFixed
}



// Report the results. 
// remember m2= m-m1


tempname cond
postfile  `cond' naive max_q sorter logor se str50 ( timepoint fix_ran uni_multi ) ///
	using ConditionalResults_`=6*`m1''mo,  replace

local j1 = 0

forval j2 = `=`m1'+1'/`m' {
	local j1 = `j1' +1

	// the contrast matrix 
	mat a = J(1, `m2', 0)
	mat a[1, `=`j2'-`m1''] = 1

	// Conditional analyses: the effect size and its variance
	mat ES  = a*condbetaFixed'
	mat VAR  = a*(condcovbetaFixed)*a'
	post `cond' (0) (`m1') (`j1') (`=el(ES, 1, 1)') ///
		(`=sqrt(el(VAR,1,1))')  ("`=`j2'*6' months | `=6*(`m1'+1)' mo") ///
		("fixed") ("multivariate")
 	
	// Naive conditional: the effect size and its variance
	mat ES  = a*naivecondbetaFixed'
	mat VAR  = a*(naivecondcovbetaFixed)*a'
	post `cond' (1) (`m1')  (`j1') (`=el(ES, 1, 1)') ///
		(`=sqrt(el(VAR,1,1))')  ("`=`j2'*6' months, moving baseline to `=6*(`m1'+1)' mo") ///
		("fixed") ("multivariate")
	
}

postclose `cond'


use ConditionalResults_`=6*`m1''mo, clear

sort naive , stable

// number of timepoints (outcomes)
// remember m2=m-m1 

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




