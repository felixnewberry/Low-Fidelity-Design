function [err_bound,err_bi, err_low] = my_airfoil_bound(Uc,Uf, N, n, r)
%Bi-fidelity error bound for diffuser

% Inputs: 
% Uc - low fidelity data n_pointsxn_samples (or course vs fine)
% Uf - high fidelity data n_pointsxn_samples
% N - total number of samples to estimate
% n - number of samples to use in bound, n << N
% r - truncation with which bi-fidelity model is found

% Outputs; 
% err_bound 
% err_bi
% err_low

rand_sample = 1:n; 

B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

% Obtain column skeleton of P
[P_s,ix] = matrixIDvR(B,r);

% Error bound inputs
normC = norm(P_s);
sb = svd(B); 
err_Bhat = norm(B-B(:,ix)*P_s);  

% % % Compute epsilon tau... 
[~, ahat_error_est,~, ~,~] = ...
    mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);

err_bound = ahat_error_est/norm(A);
err_bi = norm(A-A(:,ix)*P_s)/norm(A);
err_low = norm(Uc - Uf)/norm(Uf); 

end

