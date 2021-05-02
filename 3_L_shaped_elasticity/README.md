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

Idx_c - hmm 6 indices
Idx_f - hmm 22 indices. I suspect these are for the displacement line
x_c - course coordinates displacement line
x_f - fine coordinates displacement line
Uc, Uf, Ub, denote course, fine and bi-fidelity estimates
line - displacement line QoI
stress - stress field QoI
field - displacment field QoI

## L_design
Results of L_shaped parametric design

## fenics_inputs
xi - 800x49 
Youngs - 800x41 800 samples x 41 nodal points for course grid
inputs - python inputs

## mesh
meshes 1 through 4 

## Plots
Plot results  
