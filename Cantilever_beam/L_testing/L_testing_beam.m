clear all
close all
% clc

save_on = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;     % Line width
MS = 8;     % Marker Size
FS_leg = 16; % Font size legend


size_1 = [0,0,670,515]; 
size_2 = [0,0,1340,515]; 

FS = 28;    % Font size axis
FS_axis = 18; 
LW_axis = 1; 

% Colors
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; 
c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880]; 
c6 = [0.3010, 0.7450, 0.9330]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_samples = 100; 

Uc_nom = load('Beam_data/Uc.mat', 'Uc');
Uc_nom = Uc_nom.Uc(:,1:n_samples);
load('Beam_data/Uf.mat')
Uf = Uf(:,1:n_samples); 

Uc_opt = load('Beam_design/Uc_opt.mat');
Uc_opt = Uc_opt.Uc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply transformations + pre-processing. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = 1; 
n = r+2; 
rand_sample = 1:n; 
N = 100; 

Uc = Uc_nom; 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');
[error_bound_nom,err_Bi_nom, P_s_nom] = my_bound_bi(n,r, A, B, N);

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

% Tansform L 
Phi = A_R* pinv(B_R); 
B = Phi*B; 

[error_bound_nom2,err_Bi_nom2, P_s_nom2] = my_bound_bi(n,r, A, B, N);


Uc = Uc_opt; 
% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
[error_bound_opt,err_Bi_opt, P_s_opt] = my_bound_bi(n,r, A, B, N);

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

% Tansform L 
Phi = A_R* pinv(B_R); 
B = Phi*B; 
[error_bound_opt2,err_Bi_opt2, P_s_opt2] = my_bound_bi(n,r, A, B, N);

res_mat =  [error_bound_nom, error_bound_nom2,error_bound_opt,...
    error_bound_opt2; ...
    err_Bi_nom, err_Bi_nom2, err_Bi_opt, err_Bi_opt2]; 

res_tab = array2table(res_mat,'VariableNames',{'Nom','Nom2', 'Opt', 'Opt2'}, 'RowNames',{'Bound','Bi'})
