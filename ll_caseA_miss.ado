*! version 0.4 24mar2010
*! version 0.3 2feb2010
*! ttrikalin@mac.com  
// case A: common sigma - common correlation 
// loosely based on the programming logic in Ian White's mvmeta

prog def ll_caseA_miss

if ($be_verbose)  noi di as text _n "A: common variance - common correlation"
if ($be_verbose)  noi di as text    "             [can parse missing values]" _n

global be_verbose = 0  // do not repeat the above info at each call!

args todo b  lnf

local y $ymat
local S $Smat
local K $K
local m $m

local epsilon = 1e-6   // see below back-transforming correlations 

tempname BETA SD rho1 rho C T W dev minus2ll Wsum ll
tempname X P Wsum_miss

// get input arguments  
mat `BETA' = J(1, `m', 0)
forval i =1/`m' {
	mat `BETA'[1, `i']=el(`b',1,`i')
}
local j0 `m'

local ++j0
scalar `SD' = exp(el(`b', 1, `j0'))

local ++j0
scalar `rho1' = exp(el(`b', 1, `j0'))

// back transforming to obtain the bounded correlation 
scalar `rho' = (`m'*tanh(`rho1') + `m' -2)/(2*`m'-2)

// because of risk for numerical over- or underflow in the calculation 
// above, I am catching the correlation within espilon from its bounds. 
// i.e., -1/(m-1) and 1.  See paper
if (`rho' >=. ) {
	if (sign(`rho1')==-1) {
		scalar `rho' = -1/(`m'-1) + `epsilon'
	}
	if (sign(`rho1')==1) {
		scalar `rho' = 1 - `epsilon'
	}
}


mat `C' = J(`m', `m', `rho')
forval i=1/`m' {
	mat `C'[`i', `i']=1
}

mat `T' = `SD' * `C'

mat `Wsum' = J(`m',`m',0)
scalar `ll'= 0


forvalues i = 1/`K' {

	local hasmissing = (diag0cnt(`S'`i')>0)
	if (`hasmissing' == 0) {
		cap mat `W' = invsym(`S'`i'+`T')
		if  _rc {
			noi di as error "Problem in study `i'"
			exit -1
  		}
		mat `dev' = `y'`i'-`BETA'
		mat `minus2ll' = `m'*log(2*_pi) - log(det(`W')) + `dev' * `W' * `dev''
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
		mat  `Wsum_miss' = J(`m',`m', 0)
		forval j1=1/`nonmissing' {
			forval j2=1/`nonmissing' {
				mat `Wsum_miss'[`j1', `j2'] = `W'[`j1', `j2']
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
	scalar `ll' = `ll' - log(det(`Wsum'))/2 + `m'*log(2*_pi)/2
}
scalar `lnf' = `ll'

mat BETA = `BETA'
mat T = `T'
mat COVBETA = invsym(`Wsum')

end
