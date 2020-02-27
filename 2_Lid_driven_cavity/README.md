#Lid-driven Cavity

NavierStokesLDC.py - solve LDC problem in fenics 
run_LDC_ensemble.py - run ensemble of LDC problems
my_ldc_bound.m - calculate bound for given low and high fidelity settings
LDC_parametric_design.m - low-fidelity design process for LDC
LDC_plot.m - plot results. ie line search, response surface, error histogram

LDC_design: 
Contains results from LDC_parametric_design necessary for plotting. 
delta - the vector of pertubations to a given parameter. 
line - the corresponding error bound for a given line search
rand - random sample results to use in Polynomial Chaos Expansion. Contains the efficacy, bi-fidelity error, error bound, limits searched within and random variable realizations. 
files that contain nom and opt in them are the nominal and optimal results for a given QoI. Each one contains the bi-fidelity error, the bound, the singular values, and the bi-fidelity and low-fidelity data matrices. 

LDC_data: 
input data to for run_LDC_ensemble, saved from lDC_parametric_design which calls my_ldc_bound. 

u_meshes:
Nominal solution matrices for different QoIs and for different meshes.

Plots:
Plot of results  
