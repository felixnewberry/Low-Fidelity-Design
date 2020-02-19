function [error_bound_vec,err_Bi_vec,err_low_vec] = my_ldc_an_bound(nx, n, r, delta_Re)

1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LDC details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xi = load('LDC_data/xi_mat_2.mat');
% xi = xi.xi_2; 
% 
% entire field? - nah. 
load('u_meshes/u_64_f_2', 'u_matrix_0', 'u_matrix_1')

%%%
Uf_all = cat(3,u_matrix_0, u_matrix_1);
% This loads u_matrix 0 and 1, u mid and u vert. 


% % This loads
% % u_matrix_ 0, 1, 2, 3, 4 
% 
% % 0 is u mid v 
% % 1 is u vert u 

% % 2 is P mid
% % 3 is P vert           
% % 4 is P base

Qoi_vec = 1:2; 

error_bound_vec = zeros(length(Qoi_vec),1); 
err_Bi_vec = zeros(length(Qoi_vec),1); 
err_low_vec = zeros(length(Qoi_vec),1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate low fidelity ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delete Uc matrix file 
% delete 'u_meshes/u_matrix.mat'

% Uc_all = cat(3,Uc.u_matrix_0, Uc.u_matrix_1, Uc.u_matrix_2, Uc.u_matrix_3, Uc.u_matrix_4);
% order is u_y, u_x, p_mid, p_vert, p_base 
% Uc = Uc.u_matrix'; 

%%% Ensemble 
load('u_nu_vec_2.mat', 'nu_vec', 'u_top_vec')
% nu_vec and u_top_vec
Re_vec = u_top_vec./nu_vec;

% Apply delta

Re_vec = (1+delta_Re).*Re_vec; 

load('x_64.mat', 'x_64')
x_64 = x_64(:,1); 

x_vec = x_64; 
y_p5 = 0.5; 

1; 

Uc_all = cat(3,v_ldc(Re_vec, x_vec, y_p5), u_ldc(Re_vec, y_p5, x_vec));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fid details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% approximate N = 200 high fidelity samples using first 10. and all N low
% fid. 
% tol = 1e-4; 

% Subset of vectors for bi-fidelity error estimate
rand_sample = 1:n; 

% Step through each QoI

for i_qoi = 1:length(Qoi_vec)
    Uc = Uc_all(:,:,i_qoi)'; 
    Uf = Uf_all(:,:,i_qoi)'; 
    
    1; 
    
    % % Normalize matrices for bi-fidelity error bound
    B = Uc/norm(Uc,'fro');
    A = Uf/norm(Uf,'fro');

    B_R = B(:,rand_sample);
    A_R = A(:,rand_sample);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% L Transformation
%     Phi = A_R* pinv(B_R); 
%     B = Phi*B; 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Obtain column skeleton of P
    [P_s,ix] = matrixIDvR(B,r);

    % Error bound inputs
    normC = norm(P_s);
    sb = svd(B); 
    err_Bhat = norm(B-B(:,ix)*P_s); 
    N = 200;

    % % % Compute epsilon tau... 
    [~, ahat_error_est,~, ~,~] = ...
        mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);

    error_bound_vec(i_qoi) = ahat_error_est/norm(A);
    % efficacy = error_bound/err_Ahat;


    % Do bi-fidelity estimate to compare with bound

    err_Bi_vec(i_qoi) = norm(A-A(:,ix)*P_s)/norm(A);
%     efficacy = error_bound/err_Bi;
    err_low_vec(i_qoi) = norm(Uc - Uf)/norm(Uf); 


end

function u = u_ldc(Re,x,y)
u = Re*8*((x.^4-2*x.^3+x.^2)*(4*y.^3-2.*y))'; 
end

function v = v_ldc(Re,x,y)
v = Re.*-8*((4*x.^3-6*x.^2+2*x)*(y.^4-y.^2))'; 
end


end