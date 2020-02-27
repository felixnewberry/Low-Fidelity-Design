% Random sample 

clear all
close all 
clc

Re_nom = 100; 
Re_max = 130; 
Re_min = 70; 

u_samp = rand(200,1)*2-1;
Re_vec = Re_nom + u_samp*30;

% save('Re_vec','Re_vec')
% save('u_samp','u_samp')