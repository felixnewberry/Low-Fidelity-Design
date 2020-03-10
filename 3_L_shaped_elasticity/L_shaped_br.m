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


load('L_data/Uf_stress')
Uf = Uf(:,1:n_sim);

load('L_design/all_nom')
Uc_nom = Uc_all{3}; 
ix_nom = ix_all{3};

load('L_design/all_opt')
Uc_opt = Uc_all{3}; 
ix_opt = ix_all{3};

% random inputs
load('Fenics_inputs/xi')
% d = 7; 
% n_terms = d*d; 
xi_ref = xi(1:n_sim,:);
xi_low = xi(1:n_sim,:);

% rename and tranpose to be: Number of samples x spatial/temporal dimension of problem
Uf = Uf'; 
Uc_nom = Uc_nom'; 
Uc_opt = Uc_opt'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PC Specs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Desired polynomial order of PCE
pol = 1; 
 
% set range of number of high-fidelity samples 
N_hi = [12]; 
% N_hi = [10, 15, 20, 25];

% Set approximation rank of truncation 
r = [5]; 
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

% P is 1275 if pol = 2... for pol = 1: P=50

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

%%% calibrate with spgl1 on end point 
weights = get_matrix_weights(psi_low);
Psiweights = psi_low*weights;
% % sigma is truncation error of PCE (approximated)
sigma =  cross_val_sigma(psi_low,Uc_nom(:,10));

fprintf('Low-fidelity PC validation error of sample stress: %d %% \n', 100*sigma)

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

fprintf("r = %d \n", r);
fprintf("N_hi = %d \n", N_hi)

results_tab = array2table(100*results, 'VariableNames',{'Low','Bi'},...
    'RowNames',{'Nominal', 'Optimal'})

% P = 1, N_hi = 7, r = 4/4

% Nom = 13.7, opt = 9.86, pc = 117% 
% Nom, Opt about 8.3 for r=1... 

% 13.4 to 6.5 for r = 5, N = 7. Similar 6. 14 to 10 for 7. 
% 11.6 to 4.7 for r = 5 N = 12

% Maybe I should ramp up the number of samples for the nominal and optimal
% - ie full 800. Test effiacy on that! 
% Yes. Create new samples 

save('L_design/results_br','results', 'results_tab', 'u_bi_nom', 'u_bi_opt');
% save('L_design/results_br_all','results', 'results_tab', 'u_bi_nom', 'u_bi_opt');