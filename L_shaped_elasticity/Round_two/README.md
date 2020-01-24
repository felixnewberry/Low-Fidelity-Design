# L-shaped Elasticity 

lshape.py - L_shaped Elasticity FEniCS model
L_shaped_parametric_design - parametric low-fidelity model design
my_L_bound - Error bound calculation for given inputs
youngs_samp_gen_gaussian_cov - youngs modulus calculation
compute_eig - eigenvalue compututation used in youngs modulus
Cfunc - used in youngs modulus

## L_data
dof_coords_c.mat - Course grid x and y DOF coordinates
dof_coords_f.mat - Fine grid x and y DOF coordinates

## fenics_inputs
xi - 800x49 - why this size?

## mesh
meshes 1 through 4 

EB_Cantilever.m - beam model   
beam_parametric_design - parametric low-fidelity model design   
my_beam_bound - calclate bound for given low-fidelity model settings   
beam_plots - plot line searches, grid search, deflection realization, deflection ensemble, tip error histogram, eigenvalue decay.   

## Beam_data
Uc - course/low-fidelity vertical displacement of top chord  
Uf - fine/high-fidelity vertical displacement of top chord  
x_highfidelity - corrspondong coordinates of QoI   
xi - corresponding random inputs   

## Beam_design
Contains various design results:   
Bi_nom - nominal bi-fidelity estimate  
Bi_opt - optimal bi-fidelity estimate   
Uc_opt - optimal low-fidelity 
delta - deltas for a given search ie h1 (parameter varied) p2m1 from delta = -1 to +2 (-100 to +200 % change)
line - line search for given parameter corresponding to deltas

## Plots
Plot results  
