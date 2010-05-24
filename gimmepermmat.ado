*! version 0.2 2feb2010
*! ttrikalin@mac.com  

program define gimmepermmat , rclass

	syntax , matname(string) [ permmatname(string) ]
	
	tempname A D P X elemP
	cap mat `A' = `matname'
	if _rc noi di as error "option matname() expects a symmetric matrix. Thanks!"
	
	
	local p = colsof(`A')
	local num_full = `p' - diag0cnt(`A')
	
	// permutation matrix - start with the identity 
	mat `P' = I(`p')
	
	// if we have an empty row/column, fix the permutation matrix
	
	local goon =1 
	while(`goon') {
		local indexoffull = 0
		local indexofempty =0
		local checkforempty = 1
		local checkforfull = 1
	
		// the strategy is to construct a series of elementary permutation matrices
		mat `D' = `P'*`A'*`P''
		
		// a convenient way is to track the zeros in the diagonal 
		mat `X' = vecdiag(`D') 
	
		// find the first empty row/column and the first 
		// full row/column after the first empty one
		forval i=1/`p' {
	
			// find the first empty diagonal element 
		
			if ((el(`X',1,`i') == 0) & `checkforempty'){
				local indexofempty = `i' 
				local checkforempty = 0
			}	
			
			// only after you've found an empty search for the first full
			if ((el(`X',1,`i') != 0) & (`checkforempty'==0) & (`checkforfull'==1)) {
				local indexoffull = `i' 
				local checkforfull = 0
			}	
		}
	
		// now switch the first empty row/column with the first 
		// full (nonempty) after it
		// but only if we are still going on 
	
		if (`indexoffull' > `indexofempty') {
		
			// construct a new elementary permutation matrix
			mat `elemP'=I(`p')
			
			//move the full to the empty
			mat `elemP'[`indexofempty',`indexoffull']=1
			
			// move the empty to the full
			mat `elemP'[`indexoffull',`indexofempty']=1
			
			// get rid of the diagonal 1's (that make 
			// the row/column to "stay put")
			mat `elemP'[`indexofempty',`indexofempty']=0		
			mat `elemP'[`indexoffull',`indexoffull']=0		
			
			// get the total permutation matrix by multiplying the 
			// elementary matrices 
			mat `P' = `elemP' * `P' 
			
			
		}
		
		if (`indexoffull' ==0)  local goon=0  // we are done
	}

	if ("`permmatname'"=="")  local permmatname "P"

	// return permutation matrix	
	return matrix `permmatname' = `P'

//	mat drop `P' `D' `elemP' `A' `X'


end
