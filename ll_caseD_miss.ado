*! version 0.3 2feb2010
*! ttrikalin@mac.com  
// case D: different sigmas - autoregressive correlation 
// loosely based on the programming logic in Ian White's mvmeta

prog def ll_caseD_miss

if ($be_verbose) noi di as text _n "D: different variances - autoregressive correlation"
if ($be_verbose) noi di as text    "                         [can parse missing values]" _n

global be_verbose = 0

args todo b  lnf

local y $ymat
local S $Smat
local n $n
local p $p

local epsilon = 1e-5

tempname BETA SD rho1 rho C T W dev minus2ll Wsum ll
tempname X P Wsum_miss

// get input arguments  
mat `BETA' = J(1, `p', 0)
forval i =1/`p' {
	mat `BETA'[1, `i']=el(`b',1,`i')
}
local k `p'

mat `SD' = J(1, `p', 0)
forval i=1/`k' {
	local ++k
	mat `SD'[1, `i'] = exp(el(`b', 1, `k'))
}

local ++k
scalar `rho1' = el(`b', 1, `k')

// back transforming to obtain the bounded correlation (between 0 and 1)
scalar `rho' = 0.5*tanh(`rho1') + 0.5
if (`rho' >=. ) {
        if (sign(`rho1')==-1) {
                scalar `rho' = `epsilon'
        }
        if (sign(`rho1')==1) {
                scalar `rho' = 1 - `epsilon'
        }
}


mat `T' = diag(`SD')
forval i=1/`p' {
	forval j=1/`p' {
		if (`i' != `j') {
			mat `T'[`i', `j'] = (`rho'^abs(`i'-`j') )*sqrt(`SD'[1,`i'] * `SD'[1,`j'])
		}
	}
}


mat `Wsum' = J(`p',`p',0)
scalar `ll'= 0

forvalues i = 1/`n' {

	local hasmissing = (diag0cnt(`S'`i')>0)
	if (`hasmissing' == 0) {
		cap mat `W' = invsym(`S'`i'+`T')
		if  _rc {
			noi di as error "Problem in study `i'"
			exit -1
  		}
		mat `dev' = `y'`i'-`BETA'
		mat `minus2ll' = `p'*log(2*_pi) - log(det(`W')) + `dev' * `W' * `dev''
		mat `Wsum' = `Wsum' + `W'
	}
	if (`hasmissing'==1) {
		// get the permutation matrix to isolate the non-missing rows/columns
		gimmepermmat, matname(`S'`i')

		local nonmissing = `=colsof(`S'`i')' - diag0cnt(`S'`i')

		mat `P' = r(P)
		mat `X' = `P'*(`S'`i'+`T')*`P''
		mat `X' = `X'[1..`nonmissing', 1..`nonmissing']
		mat `W' = syminv(`X')
		mat `dev' = (`y'`i'-`BETA')*`P''
		mat `dev' = `dev'[1,1..`nonmissing']
		mat `minus2ll' = `=colsof(`S'`i')'*log(2*_pi) - log(det(`W')) + `dev' * `W' * `dev''

		// pad W with 0's for the missing rows/columns 
		mat  `Wsum_miss' = J(`p',`p', 0)
		forval j=1/`nonmissing' {
			forval k=1/`nonmissing' {
				mat `Wsum_miss'[`j', `k'] = `W'[`j', `k']
			}
		}

		// rearrange to match the row/column order of the pxp Wsum
		// and add to complete Wsum 
		mat `Wsum_miss' = `P''*`Wsum_miss'*`P'
		mat `Wsum' = `Wsum' + `Wsum_miss'
	}
	scalar `ll' = `ll' - el(`minus2ll',1,1)/2
}
if ($restricted == 1 ) {
	// of there is a study with a missing outcome, the constant part is not correct!
	scalar `ll' = `ll' - log(det(`Wsum'))/2 + `p'*log(2*_pi)/2
}

scalar `lnf' = `ll'

mat BETA = `BETA'
mat T = `T'
mat COVBETA = invsym(`Wsum')

end
