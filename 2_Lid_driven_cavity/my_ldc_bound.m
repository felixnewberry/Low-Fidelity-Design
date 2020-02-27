function [error_bound_vec,err_Bi_vec,err_low_vec] = my_ldc_bound(nx, n, r, delta_u, delta_nu_0, delta_nu_1,run_count)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LDC details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% first check standard. Have to change nu and u in python later to check
%%% new one. 
% use new RVs: twice the range

xi = load('LDC_data/xi_mat_2.mat');
xi = xi.xi_2; 

load('u_meshes/u_64_f_2', 'u_matrix_0', 'u_matrix_1', 'u_matrix_2',...
    'u_matrix_3', 'u_matrix_4')

% %%% use original random variables
% xi = load('LDC_data/xi_mat.mat');
% xi = xi.xi; 
% 
% % load high fidelity original 
% load('u_meshes/u_64_f_standard', 'u_matrix_0', 'u_matrix_1',...
%     'u_matrix_2', 'u_matrix_3', 'u_matrix_4')

%%%
Uf_all = cat(3,u_matrix_0, u_matrix_1, u_matrix_2, u_matrix_3, u_matrix_4);

% This loads
% u_matrix_ 0, 1, 2, 3, 4 

% is this correct?? Pretty vital that it is... 
% 0 is u mid v 
% 1 is u vert u 
% 2 is P mid
% 3 is P vert           
% 4 is P base

Qoi_vec = 1:5; 

error_bound_vec = zeros(length(Qoi_vec),1); 
err_Bi_vec = zeros(length(Qoi_vec),1); 
err_low_vec = zeros(length(Qoi_vec),1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate low fidelity ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save './ensemble_inputs/nx.mat' 'nx' 
save('LDC_data/inputs_vec.mat', 'nx', 'delta_u', 'delta_nu_0','delta_nu_1','run_count') ;

1; 

tic
% Call python/fenics
pyFi_ensemble = 'sudo python3.6 run_LDC_ensemble.py';
system(pyFi_ensemble); 
toc

% samps = 1:200;

% if run is succesful u_matrix is saved and bi-fidelity bound is calculated
try
% load Uc matrix
Uc= load('u_meshes/u_matrix.mat');
1; 

% delete Uc matrix file 
delete 'u_meshes/u_matrix.mat'

Uc_all = cat(3,Uc.u_matrix_0, Uc.u_matrix_1, Uc.u_matrix_2, Uc.u_matrix_3, Uc.u_matrix_4);
% order is u_y, u_x, p_mid, p_vert, p_base 
% Uc = Uc.u_matrix'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fid details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% approximate N = 200 high fidelity samples using first 10. and all N low
% fid. 
% tol = 1e-4; 



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

    % Subset of vectors for bi-fidelity error estimate 
    % Ensure column skeleton slection are used + additional samples to reach
    % total of n

    rand_sample = [1:n-length(ix),ix]; 

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
    1; 
    
    % if i_qoi == 1
%     save_label = 'u_mid';
% % elseif i_qoi == 2
% %     save_label = 'u_x';
% % elseif i_qoi == 3
% %     save_label = 'P_mid';
% % elseif i_qoi == 4
% %     save_label = 'P_vert';
% % elseif i_qoi == 5
% %     save_label = 'P_base';
% % end
% % 
% Ub = Uf(:,ix)*P_s;
% % 
% % 1; 
% % 
% save_label = 'all';

% save(strcat('LDC_design/',save_label, '_nom'),'Uc', 'Ub', 'sb')
% save(strcat('LDC_design/',save_label, '_opt'),'Uc', 'Ub', 'sb')


% save(strcat('LDC_design/',save_label, '_nom_2'),'Uc', 'Ub', 'sb')
% save(strcat('LDC_design/',save_label, '_opt_2'),'Uc', 'Ub', 'sb')
end

catch
    err_Bi_vec(:) = nan; 
    err_low_vec(:) = nan; 
end



end