function [error_bound,err_Bi, P_s] = my_bound_bi(n,r, A, B, N)
% Calculate bi-fidelity and error bound for set of low-fidelity data

rand_sample = 1:n; 
B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

% Obtain column skeleton of P
[P_s,ix] = matrixIDvR(B,r);

% Error bound inputs
normC = norm(P_s);
sb = svd(B); 
err_Bhat = norm(B-B(:,ix)*P_s); 
% N = 200;

% % % Compute epsilon tau... 
[~, ahat_error_est,~, ~,~] = ...
    mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);

error_bound = ahat_error_est/norm(A);
% efficacy = error_bound/err_Ahat;

% Do bi-fidelity estimate to compare with bound

err_Bi = norm(A-A(:,ix)*P_s)/norm(A);
%     efficacy = error_bound/err_Bi;
end

