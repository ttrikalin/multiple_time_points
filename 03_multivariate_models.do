// Perform fixed effects analyses and all REML random effects analyses 
// There are 5 variants of random effects we are fitting here with REML. They differ 
// in the way the T matrix is structured (covariance matrix of the random effect)
// case A: common sigma and common correlation (2 parameters)
// case B: different sigmas and a common correlation (5 parameters)
// case C: common sigma and a banded correlation structure (2 parameters)
// case D: different sigma and banded correlation structure (5 parameters) 
// case E: unstructruted (10 parameters) 

use data, clear
set seed 12345

local p 4

qui count 
local n = r(N)

global restricted = 1 // do REML rather than ML 


// These globals and matrices are needed to pass information to the likelihood 
// optimization programs. 

global n `n'
global ymat  "ymat"
global Smat  "Smat"
global p `p'

forval i=1/`n' {
	mat ymat`i' = J(1, `p', 0)	
	forval j =1/`p' {
		qui summ b`j' in `i' , meanonly		 
		mat ymat`i'[1,`j'] = r(mean) 
	}
	if (matmissing(ymat`i')) {
		noi di in yellow "Study `i' has missing values:"
		noi mat li ymat`i'
	}

	mat Smat`i' = J(`p', `p' , 0)
	forval k =1/`p' {
		forval m =`k'/`p' {
			qui summ  V`k'`m' in `i' , meanonly 
			mat Smat`i'[`k', `m'] = r(mean)
			mat Smat`i'[`m', `k'] = r(mean)
		}
	}	
}


// Note that study 17 has missing values. Its likelihood contribution is different
// than that of the other studies.  
// Here we replace the missing values in time point 1 and 3 (6 and 18 months)
// with zeros.  
// This will work fine for the fixed effects GLS -- see below. 
// The REML programs however will have to do some footwork to accommodate the 17th study

forval i=1/`p' {
	if (el(ymat17,1,`i')>=.)  mat ymat17[1,`i']=0
	forval j=1/`p' {
		if (el(Smat17,`i',`j')>=.)  mat Smat17[`i',`j']=0
	}
}


// All matrices must be non-singular. In this example, Smat14 (study 14) is singular
// We must assert that these matrices are positive definite; the fastest way
// would be to try-catch a Cholesky decomposition (See G Strang Linear Algebra)
// but here i will simply check the eigenvalues, as described in the paper. 
// correctmemat is a short program that does this correction quickly.

// First do the 16 studies with complete data
forval i=1/`=`n'-1' {
	correctmemat, matname(Smat`i') epsilon(0.08)
	if (r(corrected)) mat Smat`i' = r(M)
}

// now check the 17th study for the outcomes that are non-missing 
// Calculate a permutation matrix P that will permute Smat17 so that 
// the nonmissing timepoints are first and the missing follow
// this is done with a small program I wrote
gimmepermmat , matname(Smat17) 
mat P = r(P)  // this is the permutation matrix 

mat A = P*Smat17*P'  // save rearranged matrix as A
mat A= A[1..2, 1..2] // only the nonmissing values 
di det(A)	     // determinant is positive, go on
		     // if it were not, we would correct the 2x2 A matrix as above,
		     // pad the corrected 2x2 with 0's to get a 4x4
		     // and restore the order of the rows and columns 
		     // using P'*corrected_padded_matrix*P
mat drop A	     


// This is the fixed effects model, fit with GLS. The same can be fit using ML
// but I want to show this coding also. So bear with me. 
// Because i have no covariates, there is no design matrix here (it's I(4), omitted)

if (0==0) {  
	mat Wsum = J(`p', `p', 0)
        mat Wysum = J(1, `p', 0)

	// The 16 studies with all 4 time points are OK
	// turns out that for the fixed effects model,  
	// i can use the Smat17 which has 0's in the rows and columns 
	// of the missing time points; and the means vector 
	// ymat17, which also has 0's for the missing timepoints. 
	// This works out to be the same as the method by 
	// Gleser and Olkin referenced in the paper.
	// this strategy will not work OK for the REML calculations!

        forval k=1/`n' {
                mat W`k' = syminv(Smat`k')
                mat Wysum = Wysum + ymat`k' * W`k'
                mat Wsum = Wsum + W`k'
        }

        mat covbetaFixed = syminv(Wsum)
        mat betaFixed = Wysum * covbetaFixed
	
	// by definition TauFixed =0; I will use this when gathering results
	mat TauFixed = J(`p', `p', 0) 
}

// This is the fixed effects model with a different coding (REML)
// The results are identical with those from the GLS approach above

if (0==0) {
        // Fixed effects model with REML 
        program drop _all
        global be_verbose = 1  // force ll program to identify itself - avoid blunders 

        ml model d0 ll_fixed_miss (mu: b1 b2 b3 b4 , nocons), obs(`n')  collinear

        mat b0 = [betaFixed]
        ml init b0 , copy
        //ml check

        ml max , difficult iterate(30) ltolerance(1e-6)
        est store Fixed
}


// These are the five examples of random effects that are discussed in the paper
if (1==1) {
	// case A: common variance and common correlation (5 parameters)
	program drop _all 
	global be_verbose = 1  // force ll program to identify itself - avoid blunders 

	ml model d0 ll_caseA_miss (mu: b1 b2 b3 b4 , nocons) (S:) (rho:), obs(`n')  collinear 
	
	ml search S: -10 10 rho: 0 1 

	mat b0 = [betaFixed, 0, 0] 
	ml init b0 , copy 
	//ml check

	ml max , difficult iterate(30) ltolerance(1e-6)
	est store A
	mat betaA = BETA
	mat TauA = T 
	mat covbetaA = COVBETA
}

if (2==2) {
	// case B: different variances and common correlation (5 parameters)
	program drop _all 
	global be_verbose = 1  // force ll program to identify itself - avoid blunders 

	ml model d0 ll_caseB_miss (mu: b1 b2 b3 b4 , nocons) (S1:) (S2:) (S3:) (S4:) (rho:), ///
		obs(`n')  collinear 

	ml search S1: -10 10 S2: -10 10 S3: -10 10 S4: -10 10 rho: 0 1 
	
	mat b0 = [betaFixed, 1, 1,1,1,0] 
	ml init b0 , copy skip
	//ml check

	ml max , difficult iterate(30) ltolerance(1e-6)
	est store B
	mat betaB = BETA
	mat TauB = T 
	mat covbetaB = COVBETA
}

if (3==3) {
	// case C: common variance - autoregressive correlation (2 parameters)
	program drop _all 
	global be_verbose = 1  // force ll program to identify itself - avoid blunders 

	ml model d0 ll_caseC_miss (mu: b1 b2 b3 b4 , nocons) (S:) (rho:) , obs(`n')  collinear 
	
	ml search S: -10 10  rho: 0 1

	mat b0 = [betaFixed, 0, 0] 
	ml init b0 , copy skip
	//ml check

	ml max , difficult iterate(30) ltolerance(1e-6)
	est store C
	mat betaC = BETA
	mat TauC = T 
	mat covbetaC = COVBETA
}

if (4==4) {
	// case D: different variances - autoregressive correlation (5 parameters)
	program drop _all 
	global be_verbose = 1  // force ll program to identify itself - avoid blunders 

	ml model d0 ll_caseD_miss (mu: b1 b2 b3 b4 , nocons) (S1:) (S2:) (S3:) (S4:) (rho:), ///
		obs(`n')  collinear 

	ml search S1: -10 10 S2: -10 10 S3: -10 10 S4: -10 10 rho: 0 1 
	
	mat b0 = [betaFixed, 1, 1,1,1,0] 
	ml init b0 , copy skip
	//ml check

	ml max , difficult iterate(30) ltolerance(1e-6)
	est store D
	mat betaD = BETA
	mat TauD = T
	mat covbetaD = COVBETA
}

if (5==5) {
	// case E: completely unstructured (10 parameters)
	program drop _all 
	global be_verbose = 1  // force ll program to identify itself - avoid blunders 

	ml model d0 ll_caseE_miss (mu: b1 b2 b3 b4 , nocons) ///
		(S11:) (S12:) (S13:) (S14:) ///
		       (S22:) (S23:) (S24:) ///
		       	      (S33:) (S34:) ///
			             (S44:) , obs(`n')  collinear 

	ml search S11: -10 10 S12: -10 10 S13: -10 10 S14: -10 10 ///
		  	      S22: -10 10 S23: -10 10 S24: -10 10 ///
			    	          S33: -10 10 S34: -10 10 ///
				                      S44: -10 10

	
	mat b0 = [betaFixed, 1,1,1,1,1,1,1,1,1,1] 
	ml init b0 , copy skip
	//ml check

	ml max , difficult iterate(30) ltolerance(1e-6)
	est store E
	mat betaE = BETA
	mat TauE = T 
	mat covbetaE = COVBETA
}


// Use the betaX, covbetaX and TauX matrices from the multivariate models
// and reports the results. 

tempfile Results2
tempname multi
postfile  `multi' sorter logor se Q df str20 ( timepoint fix_ran uni_multi ) using `Results2', replace

local k = 0
foreach X in Fixed A B C D E {
	
	// calculate Q per model
	mat qmat = J(1, 1, 0)
	forval j=1/`n' { 
		mat qmat = qmat + (ymat`j'-beta`X')*syminv(Smat`j' + Tau`X') * ///
			(ymat`j' - beta`X')'
	}

	forval i =1/`p' {
		local k = `k' +1

		// the contrast matrix 
		mat a = J(1, `p', 0)
		mat a[1, `i'] = 1

		// the effect size and its variance
		mat ES  = a*beta`X''
		mat VAR  = a*(covbeta`X')*a'

		post `multi' (`k') (`=el(ES, 1, 1)') ///
		(`=sqrt(el(VAR,1,1))') (`=el(qmat,1,1)') ///
			(`=`p'*(`n'-1)') ("`=`i'*6' months") ("`X'") ("multivariate")
		
	}
}

postclose `multi'


use Results, clear // we saved the univariate results here
append using `Results2'
save Results, replace 
