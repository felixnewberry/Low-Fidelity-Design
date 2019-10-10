function [error_bound,err_Bi,efficacy] = my_ldc_bound(QoI,nx, n, r, delta_u, delta_nu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LDC details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi = load('LDC_data/xi_mat.mat');
xi = xi.xi_mat; 

% 0 u_mat, 1 P_mat, 2 u_vec mid, 3 p_vec mid, 4 p_vec top

% u_matrix_32 is 66 (x and y? should just be 33... p_field_matrix is 1089 
% p_line_matrix is 33

% QoI 2 works. 
if QoI == 0
    Uf= load('u_meshes/u_field_matrix_32.mat');
    Uf = Uf.u_field_matrix';
elseif QoI == 1
    Uf= load('u_meshes/p_field_matrix_32.mat');
    Uf = Uf.p_field_matrix_32';
elseif QoI == 2
    Uf= load('u_meshes/u_matrix_32.mat');
    Uf = Uf.u_matrix';
elseif QoI == 3
    Uf= load('u_meshes/p_line_matrix_mid_32.mat');
    Uf = Uf.p_line_matrix';
elseif QoI == 4
    % this is not right.
    Uf= load('u_meshes/p_line_matrix_top_32.mat');
    Uf = Uf.p_line_matrix';
end 
    
% % load('home/felixnewberry/Documents/Research/10_low_fidelity_design/ensemble_inputs/x_all.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate low fidelity ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save './ensemble_inputs/nx.mat' 'nx' 
save('LDC_data/inputs_vec.mat', 'nx', 'delta_u', 'delta_nu','QoI') ;

1; 

tic
% Call python/fenics
pyFi_ensemble = 'sudo python3.6 run_LDC_ensemble.py';
system(pyFi_ensemble); 
toc

% samps = 1:200;

Uc= load('u_meshes/u_matrix.mat');
Uc = Uc.u_matrix'; 

1; 

% % save("Uc_better","Uc");
% 
% % Uf = Uf(:,samps);
% % Uc = Uc(:,samps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fid details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% approximate N = 200 high fidelity samples using first 10. and all N low
% fid. 
% tol = 1e-4; 

% Subset of vectors for bi-fidelity error estimate
rand_sample = 1:n; 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);


% are these steps correct???? Hmm... 
% Possibly the inputs are calculated using all low-fidelity data... I think
% so... 

% This means I should git push and re run the beam scripts... -rank 1 so
% not likely to change much but check all the same.

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

error_bound = ahat_error_est/norm(A);
% efficacy = error_bound/err_Ahat;


% Do bi-fidelity estimate to compare with bound

err_Bi = norm(A-A(:,ix)*P_s)/norm(A);
efficacy = error_bound/err_Bi;

if QoI == 0
    save_label = 'u_field';
elseif QoI == 1
    save_label = 'P_field';
elseif QoI == 2
    save_label = 'u_mid';
elseif QoI == 3
    save_label = 'P_mid';
elseif QoI == 4
    save_label = 'P_top';
end

Ub = Uf(:,ix)*P_s;
save(strcat('LDC_design/',save_label, '_nom'),'Uc', 'Ub', 'sb', 'error_bound','err_Bi')

% save(strcat(save_label, '_opt',Uc, Ub, sb, error_bound,err_Bi)

% % save('Uc_opt','Uc')
% % save('Ub_opt','Ub')
% % save('sb_opt','sb')

% % 
% % % debug plot
% % err_Bhat
% % % 9.7e-4 for case 2 (good)
% % % 2.8e-2 for case 1 (bad)
% % error_bound
% % 
% % figure
% % plot(Uc)
% % 
% % figure
% % semilogy(sb)
% % % rank 6 ... I can probably reduce R alot. 
% % 
% % 
% % 
% % 1; 
% 
% % Ub = Uf(:,ix)*P_s;
% % 
% % save('Uc_opt','Uc')
% % save('Ub_opt','Ub')
% % save('sb_opt','sb')
% % save('bound_opt','error_bound')
% % save('Ahat_opt','err_Ahat')
% 
% % save('Uc_nom','Uc')
% % save('Ub_nom','Ub')
% % save('sb_nom','sb')
% % save('bound_nom','error_bound')
% % save('Ahat_nom','err_Ahat')
% 
% % save a bunch of stuff, ie Uc, Bi, Sb _nom and _opt
% 
1; 


end