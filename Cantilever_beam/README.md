# Cantilever Beam 

EB_Camtilever.m - beam model   
epsilon_tau_parameter - parametric low-fidelity model design   
my_beam_bound - calclate bound for given low-fidelity model settings   
beam_plots - plot deflection realization and ensemble, line searches  

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
