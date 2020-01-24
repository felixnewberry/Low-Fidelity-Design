clear all 
close all
clc

%%% Test the error bound, why is it failing?


nsim = 200; 

r = 6; 
n = r+10; % subsample A and B


rand_sample = 1:n; 


% Stress
load('Uf_stress')
load('Uc_stress')

samps = 1:nsim; 

Uf = Uf(:,samps);
Uc = Uc(:,samps);

% % Normalize matrices
% Uc= Uc;
% Uf = Uf;

% Bi-fid matrices. 
% I think this normalization step is broken, I normalize by the entire
% dataset, not the subset... See LDC for reference. 

% B_R = Uc(:,rand_sample)/norm(Uc,'fro');
% A_R = Uf(:,rand_sample)/norm(Uf,'fro');

% Like this: 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

% Obtain column skeleton of P
[P_s,ix] = matrixIDvR(B,n);

% Inputs
normC = norm(P_s);
sb = svd(B); 
err_Bhat = norm(B-B(:,ix)*P_s); 
N = nsim;

% % % Compute epsilon tau... 
[~, ahat_error_est,~, ~,~] = ...
    mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);

error_bound = ahat_error_est/norm(A);

err_Bi = norm(A-A(:,ix)*P_s)/norm(A);
%     efficacy = error_bound/err_Bi;
% err_low = norm(Uc - Uf)/norm(Uf); 

error_bound
err_Bi

% Mistake is fixed!! Now run samples anew! First - go to dog park :) 