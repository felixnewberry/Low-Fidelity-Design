% Basis reduction test: 

% Question: Does an alternative bi-fidelity method work better with the
% optimal low-fidelity model that was found? 

clear all
close all
clc

%%% To do: 

% Test case works well. 
% Improve code: Try out for different N_hi and p and save results. 

% Write now I don't have a measure of the efficey of p choice. 
% I should report some polynomial order for tuning. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tip error, or entire deflection? 

n_sim = 200; 

load('u_meshes/u_64_f_2', 'u_matrix_0', 'u_matrix_1', 'u_matrix_2',...
    'u_matrix_3', 'u_matrix_4')

% % U mid  -fails
% Uf = u_matrix_0'; 
% load('LDC_design/Nom_u_mid')
% Uc_nom = Uc;
% ix_nom = rand_sample(1:r);
% load('LDC_design/Opt_u_mid')
% Uc_opt = Uc;
% ix_opt = rand_sample(1:r); 

% P base - works (with worse pce error... 
Uf = u_matrix_4'; 
load('LDC_design/Nom_p_base')
Uc_nom = Uc;
ix_nom = rand_sample(1:r);
load('LDC_design/Opt_p_base')
Uc_opt = Uc;
ix_opt = rand_sample(1:r); 


% random inputs
load('LDC_data/xi_mat_2')
xi = xi_2; 
% xi_ref dim: Number of reference samples x stochastic dimension d
xi_ref = xi(1:n_sim,:); 
xi_low = xi_ref; 


% rename and tranpose to be: Number of samples x spatial/temporal dimension of problem
Uf = Uf'; 
Uc_nom = Uc_nom'; 
Uc_opt = Uc_opt'; 

% Adjust data: no point trying to approximate 0s at edges. 
Uf = Uf(:,2:end-1); 
Uc_nom = Uc_nom(:,2:end-1); 
Uc_opt = Uc_opt(:,2:end-1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PC Specs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Desired polynomial order of PCE
pol = 3; 
 
% set range of number of high-fidelity samples 
N_hi = [2]; 
% N_hi = [10, 15, 20, 25];

% Set approximation rank of truncation 
r = [1]; 
% The covariance matrix is always rank 1 - anything else will break for
% beam. 

% tolerance on residual used in spgl1 solver
sigma = .001;

index_pc = nD_polynomial_array(size(xi_ref,2),pol); 
% assembles index of polynomials for a given stochastic dimension and 
% polynomail order. 
% one column for each variable (ie dimension d). Number of rows for set of
% P basis functions based on order pol (or p) ie 0 0, 1 0, 0 1 for p = 1
% yielding P = 3

% size of PC basis (set of P basis functions 
P = size(index_pc,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Hi-fid PC basis - in practice not available, just a limited number. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% measurement matrix used by the reference model
% size: number of data points x number of basis functions

psi_ref = zeros(size(xi_ref,1),P);

for isim=1:size(xi_ref,1)
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_ref = piset(xi_ref(isim,:),index_pc);
    psi_ref(isim,:) = crow_ref(1:P);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Low-fid Solution - nominal and optimal. - tune p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% measurement matrix used by low-fidelity model
psi_low = zeros(size(xi_low,1),P);

for isim=1:size(xi_low,1)
    % piset evaluates legendre polynomials at xi.
    crow_low = piset(xi_low(isim,:),index_pc);
    psi_low(isim,:) = crow_low(1:P);
end

opts = spgSetParms('iterations',10000,'verbosity',0);


% low fidelity coefficients solved via l_1 minimization
c_low_nom = spg_mmv(psi_low, Uc_nom, sigma*norm(Uc_nom),opts);
% transposed for convenience
c_low_nom = c_low_nom';

% low fidelity coefficients solved via l_1 minimization
c_low_opt = spg_mmv(psi_low, Uc_opt, sigma*norm(Uc_opt),opts);
% transposed for convenience
c_low_opt = c_low_opt';

% Tune p based of low.. 
% Alternative way to arrive at c_low: 
opts = spgSetParms('iterations',10000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);

%%% calibrate with spgl1 for one point - point 20 
weights = get_matrix_weights(psi_low);
Psiweights = psi_low*weights;
% % sigma is truncation error of PCE (approximated)

sigma_val =  cross_val_sigma(psi_low,Uc_nom(:,20));

fprintf('Low-fidelity PC validation error: %d %% \n', 100*sigma_val)
1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make elegant and test for different N. Whoop! This is a result. 
        
[error_low_nom, error_bi_nom, u_bi_nom] = my_br_estimate(r, N_hi, sigma, Uc_nom, Uf, psi_ref, c_low_nom, ix_nom);
[error_low_opt, error_bi_opt, u_bi_opt] = my_br_estimate(r, N_hi, sigma,  Uc_opt, Uf, psi_ref, c_low_opt, ix_opt);

%%% Print table of results
fprintf('Results basis reduction: \n');
results = [error_low_nom, error_bi_nom; ...
    error_low_opt, error_bi_opt]; 

results_tab = array2table(100*results, 'VariableNames',{'Low','Bi'},...
    'RowNames',{'Nominal', 'Optimal'})


% save('LDC_design/results_br','results', 'results_tab', 'u_bi_nom', 'u_bi_opt');