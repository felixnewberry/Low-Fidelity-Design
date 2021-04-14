function [error_bound_vec,err_Bi_vec,err_low_vec] = my_ldc_bound(nx, n, r, delta_u, delta_nu, nom_opt)

%%% Inputs

% Revisit inputs... 

% X - vector of deltas 
% X(1) = delta_u;
% X(2) = delta_nu;
% nx - grid resultion of FEniCS code 
% n - truncation - subset of data from which bound is estimated. 
% r - truncation with which bi-fidelity model is found

% mode:
% 0 test u and nu
% 1 save nominal and optimal results

%%% Outputs
% error_bound 
% err_Ahat 
% efficiacy 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LDC details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% save('LDC_data/inputs_vec.mat', 'nx', 'delta_u', 'delta_nu_0','delta_nu_1','run_count') ;
% save('LDC_data/inputs_vec.mat', 'nx', 'delta_u', 'delta_nu_0','delta_nu_1') ;
% delta_u = X(1); 
% delta_nu = X(2); 

1; 

save('LDC_data/inputs_vec.mat', 'nx', 'delta_u', 'delta_nu');

tic
% Call python/fenics
% pyFi_ensemble = 'sudo python3.6 run_LDC_ensemble.py';
pyFi_ensemble = 'sudo python3 run_LDC_ensemble.py';

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

    % Obtain column skeleton of P
    [P_s,ix] = matrixIDvR(B,r);
    
    % Error bound inputs
    normC = norm(P_s);
    sb = svd(B); 
    err_Bhat = norm(B-B(:,ix)*P_s); 
    N = 200;
    
    % Subset of vectors for bi-fidelity error estimate 
    % Ensure column skeleton slection are used + additional samples to reach
    % total of n
    rand_sample = [ix, getfield(setxor(ix,1:N), {1:n-numel(ix)})];

    B_R = B(:,rand_sample);
    A_R = A(:,rand_sample);

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
    
    
    if nom_opt == 1
        if i_qoi == 1
            % U Mid - first qoi 
            Ub = Uf(:,ix)*P_s;
            save('LDC_design/Nom_u_mid','Uc', 'Ub', 'sb','rand_sample','n','r')
        elseif i_qoi == 5
            Ub = Uf(:,ix)*P_s;
            save('LDC_design/Nom_p_base','Uc', 'Ub', 'sb','rand_sample','n','r')
        end
    elseif nom_opt == 2
        if i_qoi == 1
            Ub = Uf(:,ix)*P_s;
            save('LDC_design/Opt_u_mid','Uc', 'Ub', 'sb','rand_sample','n','r')
        end
    elseif nom_opt == 3
        if i_qoi == 5
            Ub = Uf(:,ix)*P_s;
            save('LDC_design/Opt_p_base','Uc', 'Ub', 'sb','rand_sample','n','r')    
        end
    end
    
end

if nom_opt == 1
    save('LDC_design/Nom_errors_all','error_bound_vec','err_Bi_vec','err_low_vec','n','r')
end

catch
    err_Bi_vec(:) = nan; 
    err_low_vec(:) = nan; 
    1;
    
end

end