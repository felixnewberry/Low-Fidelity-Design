function [error_bound_all,error_Bi_all] = my_L_bound(X,nsim, n, r,...
    mode_delta,n_reps)


%%% Inputs

% X - vector of deltas 
% nsim - number of runs within one sample, ie 100 for cantilever beam
% n - truncation - subset of data from which bound is estimated. 
% r - truncation with which bi-fidelity model is found

% mode delta:
% 0 is nu
% 1 is E
% 2 is corr length
% 3 is sigma
% 4 is theta - changing trapezoid angle of slope of trapezoid. Nominal value 0. Max
% maybe pi/8 to start. (corresponds to 0.4142 in max change to applied load)
% 5 is q: changing mean width of trapezoid
% 6 nu, corr and sigma

% mode_qoi: 
% 0 is line from x = 0, y = 0 to y = 1 (change r to 6)
% 1 is displacement field 
% 2 is stress field 

% compute all three qoi each time and just list results for all. 

% n_reps number of repetitions of bound samples

%%% Outputs
% error_bound - average over repetitions
% err_Ahat - average over repetitions
% efficiacy - average over repetionns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% L-shaped elaticity details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input: Define the grid levels
coarse_level = 1; % 1-4my_L_bound

% Number of samples
% nsim = 800;
% randn('state',1) %xi already generated. 

% Define sigma, correlation length, and dimension of the one-D Gaussian defining youngs

% % Test 1 
% sigma = .85;
% corr_length = 0.1;

% Test 2 
sigma = .85;
corr_length = 0.1;

d = 7;
youngs_0 = 0.1; % was 0.1. not sure what to make of this. 
num_kl_grid = 8 * d; % Note this should be appropriately checked

nu = 0.3; 
delta_y = 0;
delta_q = 0; 
theta = 0; 

if mode_delta == 0
    nu = nu*(1+X); 
elseif mode_delta ==1
    delta_y = X; 
%     delta_y = 0; 
elseif mode_delta == 2
    corr_length = corr_length*(1+X); 
elseif mode_delta == 3
    sigma = sigma*(1+X); 
elseif mode_delta == 4
    % limit should be pm pi/8 to start (maybe less)
    theta = X; 
elseif mode_delta == 5
    delta_q = X; 
elseif mode_delta == 6
    nu = nu*(1+X(1)); 
    corr_length = corr_length*(1+X(2)); 
    sigma = sigma*(1+X(3)); 
    theta = X(4); 
end
  
1; 

save('./fenics_inputs/inputs.mat','nu','theta','delta_q')

% Extract number of degrees of freedom and nodal coordinates
frid = fopen(['./mesh/mesh_' num2str(coarse_level) '.txt'],'r');
mesh_info = fscanf(frid,'%d',[1,2]);
xy_coarse = fscanf(frid,'%g %g',[2 mesh_info(1)]); xy_coarse = xy_coarse';

% Set up course grid for each qoi. 
coarse_grid_0 = 6; 
coarse_grid_1 = size(xy_coarse,1);
coarse_grid_2 = size(xy_coarse,1);


% % xi 
load('fenics_inputs/xi')

%%%%%%%%%%%%COARSE GRID
%2A generate youngs data
youngs_samp_gen_gaussian_cov(youngs_0, sigma, corr_length, d, xy_coarse, xi, num_kl_grid,delta_y)

% Generate mesh file name
system(['cp ./mesh/mesh_' num2str(coarse_level) '.xml' ' ./mesh/mesh.xml'])

%3A run python
tic;
system('sudo python3.6 lshape.py')
t_c = toc/nsim;
1; 

%4A upload data to matrix
Uc_0=zeros(coarse_grid_0,nsim);    
Uc_1=zeros(coarse_grid_1,nsim);    
Uc_2=zeros(coarse_grid_2,nsim);    


for i=1:nsim
    filename = [ 'L_data/solution_0.' num2str(i-1) '.txt' ];
    fileID = fopen(filename,'r');
    Uc_0(:,i)=fscanf(fileID, '%f');
    fclose(fileID);
    
    filename = [ 'L_data/solution_1.' num2str(i-1) '.txt' ];
    fileID = fopen(filename,'r');
    Uc_1(:,i)=fscanf(fileID, '%f');
    fclose(fileID);
    
    filename = [ 'L_data/solution_2.' num2str(i-1) '.txt' ];
    fileID = fopen(filename,'r');
    Uc_2(:,i)=fscanf(fileID, '%f');
    fclose(fileID);
end

% % save course data
% file = sprintf( 'L_data/U%d.mat',coarse_level);
% save(file,'U','-mat')

system('rm -f L_data/solution_*');

%Uf fine

% indices for identifying line from x = 0, y = 0 to y = 1. 
load('L_data/Idx_f', 'Idx_f')
load('L_data/Idx_c', 'Idx_c')
load('L_data/Uf_line', 'U')
Uf_0 = U(Idx_f,:); 

Uc_0 = Uc_0(Idx_c,:); 

load('L_data/Uf_field', 'U')
Uf_1 = U; 


load('L_data/Uf_stress', 'Uf')
Uf_2 = Uf; 


% Group into Uf_all and Uc_all

Uc_all = {Uc_0, Uc_1, Uc_2}; 
Uf_all = {Uf_0, Uf_1, Uf_2}; 

samps = 1:nsim; 

error_bound_all = zeros(3,1); 
error_Bi_all = zeros(3,1); 

Ub_all = cell(1,3);
sb_all = cell(1,3);

for i_qoi = 1:3
    Uf = Uf_all{i_qoi}; 
    Uc = Uc_all{i_qoi}; 
    
    Uf = Uf(:,samps);
    Uc = Uc(:,samps);

%     1; 
%     mean(Uc(:))


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Bi-fid details
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     error_bound_vec = zeros(1,n_reps); 

%     for i_reps = 1:n_reps

    % Subset of vectors for bi-fidelity error estimate
    % rand_sample = randsample(nsim,n);
    rand_sample = 1:n; 

    % % Normalize matrices for bi-fidelity error bound
    B = Uc/norm(Uc,'fro');
    A = Uf/norm(Uf,'fro');

    B_R = B(:,rand_sample);
    A_R = A(:,rand_sample);

    % Obtain column skeleton of P
    [P_s,ix] = matrixIDvR(B,r);

    % Inputs
    normC = norm(P_s);
    sb = svd(B); 
    err_Bhat = norm(B-B(:,ix)*P_s); 
    N = nsim;

    % % % Compute epsilon tau... 
    [~, ahat_error_est,~, ~,~] = ...
        mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);

    error_bound_all(i_qoi) = ahat_error_est/norm(A);

%     end

    error_Bi_all(i_qoi) = norm(A-A(:,ix)*P_s)/norm(A);
    
    
    Ub_all{i_qoi} = Uf(:,ix)*P_s;
    sb_all{i_qoi} = sb; 

    
end

save_label = 'all';
% 
save(strcat('L_design/',save_label, '_nom'),'Uc_all', 'Ub_all', 'sb_all', 'error_bound_all', 'error_Bi_all');

% save(strcat('L_design/',save_label, '_opt'),'Uc_all', 'Ub_all', 'sb_all', 'error_bound_all', 'error_Bi_all')

end