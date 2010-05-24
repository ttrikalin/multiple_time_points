// Stata code for main analyses described in Trikalinos and Olkin
// Meta-analysis of outcomes at multiple follow up times: a multivariate approach
// 

// First enter the data into the software 
// nAB is the numerator for arm A and time point B
// NA is the denominator for arm A

// the 17th study has no data for time points 1 and 3

clear 
input n11 n12 n13 n14 N1 n21 n22 n23 n24 N2 
//n11	n12	n13	n14	N1	n21	n22	n23	n24	N2	 
16	11	4	4	19	20	12	8	3	22	 
22	18	15	15	34	22	12	8	6	35	 
44	21	10	3	72	40	15	3	0	68	 
19	14	5	2	22	12	5	4	3	20	 
62	42	26	15	70	27	13	6	5	32	 
130	80	47	30	183	65	33	14	11	94	 
24	13	5	3	26	30	18	10	9	50	 
51	37	19	11	61	44	30	19	15	55	 
30	23	13	10	36	17	12	4	4	25	 
43	19	8	6	45	35	14	4	0	35	 
169	106	67	51	246	139	76	42	35	208	 
279	170	97	73	386	97	46	21	8	141	 
56	34	21	20	59	30	17	9	7	32	 
42	18	9	9	45	10	3	1	1	15	 
14	13	12	9	14	18	14	13	12	18	 
21	12	6	5	26	15	10	4	1	19	 
.	42	.	23	74	.	40	.	30	75	 
end 


// effect sizes (log odds ratios) bi; i-th time point
gen b1 = logit(n11/N1) - logit(n21/N2)  
gen b2 = logit(n12/N1) - logit(n22/N2)  
gen b3 = logit(n13/N1) - logit(n23/N2)  
gen b4 = logit(n14/N1) - logit(n24/N2)  

// make a continuity correction by adding 0.5 to all cells 
// for calculating b1 and b4 in studies 10 and 15
replace b1 = logit((n11+.5)/(N1+1)) - logit((n21+.5)/(N2+1)) if b1 ==.   
replace b4 = logit((n14+.5)/(N1+1)) - logit((n24+.5)/(N2+1)) if b4 ==.   

// variances Vii; i-th time point
gen V11 = 1/(n11) + 1/(N1-n11) + 1/n21 + 1/(N2-n21) 
gen V22 = 1/(n12) + 1/(N1-n12) + 1/n22 + 1/(N2-n22) 
gen V33 = 1/(n13) + 1/(N1-n13) + 1/n23 + 1/(N2-n23) 
gen V44 = 1/(n14) + 1/(N1-n14) + 1/n24 + 1/(N2-n24) 

// continuity correction for V11 and V44 in studies 10 and 15
replace V11 = 1/(n11+.5) + 1/(N1-n11+.5) + 1/(n21+.5) + 1/(N2-n21+.5) if V11==.
replace V44 = 1/(n14+.5) + 1/(N1-n14+.5) + 1/(n24+.5) + 1/(N2-n24+.5) if V44==. 

// covariances Vij; between i-th and j-th time point
gen V12 = N1/(n11 * (N1 - n12)) + N2/(n21 * (N2-n22)) 
gen V13 = N1/(n11 * (N1 - n13)) + N2/(n21 * (N2-n23))  
gen V14 = N1/(n11 * (N1 - n14)) + N2/(n21 * (N2-n24)) 
gen V23 = N1/(n12 * (N1 - n13)) + N2/(n22 * (N2-n23)) 
gen V24 = N1/(n12 * (N1 - n14)) + N2/(n22 * (N2-n24)) 
gen V34 = N1/(n13 * (N1 - n14)) + N2/(n23 * (N2-n24)) 

// continuity correction for V12 in study 15
replace V12 = (N1+1)/( (n11+.5) * (N1 - n12+.5) ) + (N2+1)/((n21+.5) * (N2-n22+.5)) if V12==. 

save data, replace 

