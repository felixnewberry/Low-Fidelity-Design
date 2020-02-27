clear all
close all
clc
% bound test 

load('L_data/Uc_stress.mat')
load('L_data/Uf_stress.mat')

r = 6; 
n = r+10; 
nsim = n; 

rand_sample = 1:nsim; 

B_R = Uc(:,rand_sample)/norm(Uc,'fro');
A_R = Uf(:,rand_sample)/norm(Uf,'fro');

[P_s,ix] = matrixIDvR(B_R,n);

% Inputs
normC = norm(P_s);
sb = svd(B_R); 
err_Bhat = norm(B_R-B_R(:,ix)*P_s); 
N = nsim;

% % % Compute epsilon tau... 
[delta_eps, ahat_error_est,Ir, min_de1,min_de2] = ...
    mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);

ahat_error_est/norm(A_R)

[P_s_r,ix_r] = matrixIDvR(Uc,r);
err_Ahat = norm(Uf-Uf(:,ix_r)*P_s_r)/norm(Uf)
% efficacy = error_bound/err_Ahat;
