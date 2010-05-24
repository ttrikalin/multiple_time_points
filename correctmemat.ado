*! version 0.3 1mar2010
*! version 0.2 2feb2010
*! ttrikalin@mac.com  

program define correctmemat , rclass
version 8.2 

syntax , matname(string) [ epsilon(numlist >1e-8 <2) ]

	tempname  x v M
	if ("`epsilon'"=="") local epsilon = 1e-5 

        cap mat `M' = `matname'
        if _rc noi di as error "option matname() expects a symmetric matrix. Thanks!"

	mat symeigen `x' `v' = `M'
        local invertible = 1
        forval j=1/`=colsof(`M')' {
                local lamda = el(`v', 1, `j')
                if (`lamda' <= `epsilon') {
                        local invertible = 0
                }
        }
        if (`invertible'==0) {
		mat `v'= `v' + J(1, `=colsof(`M')', `epsilon')
                mat `M' = `x'*diag(`v')*`x''
                noi di as text _n "Matrix `matname' was (near) singular and is now corrected" _n
        }

	return local epsilon = `epsilon'
	return scalar corrected = `=1-`invertible''
	return matrix M = `M'
end
