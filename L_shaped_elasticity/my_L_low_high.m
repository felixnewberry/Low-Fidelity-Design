function [error_bound,err_Ahat, Uc, Uf] = my_L_low_high(X,nsim, n, r,n_reps,...
    mode_delta, mode_qoi)

%%% Inputs

% % % X - vector of deltas 
% % % nsim - number of runs within one sample, ie 100 for cantilever beam
% % % n - truncation - subset of data from which bound is estimated. 
% % % r - truncation with which bi-fidelity model is found
% % 
% % % mode:
% % % 0 test w
% % % 1 test t1
% % % 2 test t2 
% % % 3 test t3
% % % 4 test 
% % 
% % % n_reps number of repetitions of bound samples
% % 
% % %%% Outputs
% % % error_bound - average over repetitions
% % % err_Ahat - average over repetitions
% % % efficiacy - average over repetionns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% L-shaped elaticity details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input: Define the grid levels
coarse_level = 1; % 1-4
fine_level = 4; 

% Number of samples
% nsim = 800;
% randn('state',1) %xi already generated. 

% Define sigma, correlation length, and dimension of the one-D Gaussian defining youngs
sigma = .85;
corr_length = 0.1;
d = 7;
youngs_0 = 0.1; % was 0.1. not sure what to make of this. 
num_kl_grid = 8 * d; % Note this should be appropriately checked

nu = 0.3; 
delta_y = 0; 

if mode_delta == 0
    nu = nu*(1+X); 
elseif mode_delta ==1
    delta_y = X; 
end
save('./fenics_inputs/inputs.mat','nu','mode_qoi')

% Extract number of degrees of freedom and nodal coordinates
frid = fopen(['./mesh/mesh_' num2str(coarse_level) '.txt'],'r');
mesh_info = fscanf(frid,'%d',[1,2]);
xy_coarse = fscanf(frid,'%g %g',[2 mesh_info(1)]); xy_coarse = xy_coarse';
% full field
if mode_qoi == 0
    coarse_grid = 6; 
elseif mode_qoi == 1
    coarse_grid = size(xy_coarse,1);
end


frid = fopen(['./mesh/mesh_' num2str(fine_level) '.txt'],'r');
mesh_info = fscanf(frid,'%d',[1,2]);
xy_fine = fscanf(frid,'%g %g',[2 mesh_info(1)]); xy_fine = xy_fine';
if mode_qoi == 0
    fine_grid = 22; 
elseif mode_qoi == 1
    fine_grid = size(xy_fine,1);
end

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
Uc=zeros(coarse_grid,nsim);    
for i=1:nsim
    filename = [ 'L_data/solution.' num2str(i-1) '.txt' ];
    fileID = fopen(filename,'r');
    Uc(:,i)=fscanf(fileID, '%f');
    fclose(fileID);
end

% % save course data
% file = sprintf( 'L_data/U%d.mat',coarse_level);
% save(file,'U','-mat')

system('rm -f L_data/solution.*');

% %Uf fine
% load('L_data/Uf')
% Uf = U; 

youngs_samp_gen_gaussian_cov(youngs_0, sigma, corr_length, d, xy_fine, xi, num_kl_grid,delta_y)

% Generate mesh file name
system(['cp ./mesh/mesh_' num2str(fine_level) '.xml' ' ./mesh/mesh.xml'])

%3A run python
tic;
system('sudo python3.6 lshape.py')
t_c = toc/nsim;
1; 

Uf=zeros(fine_grid,nsim);    
for i=1:nsim
    filename = [ 'L_data/solution.' num2str(i-1) '.txt' ];
    fileID = fopen(filename,'r');
    Uf(:,i)=fscanf(fileID, '%f');
    fclose(fileID);
end

% % save course data
% file = sprintf( 'L_data/U%d.mat',coarse_level);
% save(file,'U','-mat')

system('rm -f L_data/solution.*');


samps = 1:nsim; 

Uf = Uf(:,samps);
Uc = Uc(:,samps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fid details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_bound_vec = zeros(1,n_reps); 

for i_reps = 1:n_reps

% Subset of vectors for bi-fidelity error estimate
% rand_sample = randsample(nsim,n);
rand_sample = 1:n; 

% % Normalize matrices
% Uc= Uc;
% Uf = Uf;

% Bi-fid matrices. 
B_R = Uc(:,rand_sample)/norm(Uc,'fro');
A_R = Uf(:,rand_sample)/norm(Uf,'fro');

% % % Normalize matrices
% B_R= B_R/norm(B_R,'fro');
% A_R = A_R/norm(A_R,'fro');

% Obtain column skeleton of P
% In the paper they use 100 samples of Uc here... 
[P_s,ix] = matrixIDvR(B_R,n);

% Inputs
normC = norm(P_s);
sb = svd(B_R); 
err_Bhat = norm(B_R-B_R(:,ix)*P_s); 
N = nsim;


% % % Compute epsilon tau... 
[delta_eps, ahat_error_est,Ir, min_de1,min_de2] = ...
    mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);

error_bound_vec(i_reps) = ahat_error_est/norm(A_R);

end

error_bound = mean(error_bound_vec); 

[P_s_r,ix_r] = matrixIDvR(Uc,r);
err_Ahat = norm(Uf-Uf(:,ix_r)*P_s_r)/norm(Uf);
efficacy = error_bound/err_Ahat;


% tip error calculation
errors_Ahat = (Uf-Uf(:,ix_r)*P_s_r); %/norm(Uf);
% tip_error_1 = errors_Ahat(end,:)/norm(Uf(end,:)); 

% or 
% tip_error = errors_Ahat(end,:)./Uf(end,:); 

% save('tip_error_regular','tip_error')
% save('tip_error_optimized','tip_error')

% Bi = Uf(:,ix_r)*P_s_r; 
% % save('Beam_design/Bi_nom','Bi')
% save('Beam_design/Bi_opt','Bi')

1;


end