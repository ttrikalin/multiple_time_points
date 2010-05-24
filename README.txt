// Stata code for main analyses described in Trikalinos and Olkin
// Meta-analysis of outcomes at multiple follow up times: 
// a multivariate approach
// 
// Please save all .ado files in Stata's  working directory 
// (or anywhere in the path)
//
// Then run 
// - 01_prepare_data.do to prepare the data
// - 02_univariate_models.do to save univariate analyses (you'll need metan)
// - 03_multivariate_models.do for fixed and random multivariate analyses 
//   (you'll need all the .ado files with 
//     -likelihood programmes: ll_caseX_miss.ado, [X="fixed",or "A" to "E"] 
//     -correctmemat.ado: "corrects" singular matrices to a near posdef
//     -gimmepermmat.ado: returns a permutation matrix that I use to 
//      isolate non-missing time-points -- see code)
// - 04a_perform_conditional_analyses_fixed.do for an example of conditional analyses on 6 months with fixed effects  
//
// Thanks, 
//
// tom trikalinos
// ttrikalin@mac.com | ttrikalinos@tuftsmedicalcenter.org 
// 
