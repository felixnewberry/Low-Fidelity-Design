function [error_bound,err_Ahat,efficacy] = my_L_bound(X,nsim, n, r,...
    mode_delta,n_reps,mode_qoi)


%%% Inputs

% X - vector of deltas 
% nsim - number of runs within one sample, ie 100 for cantilever beam
% n - truncation - subset of data from which bound is estimated. 
% r - truncation with which bi-fidelity model is found

% mode:
% 0 is nu
% 1 is E
% 2 is corr length
% 3 is sigma
% 4 is theta - changing trapezoid angle of slope of trapezoid. Nominal value 0. Max
% maybe pi/8 to start. (corresponds to 0.4142 in max change to applied load)
% 5 is q: changing mean width of trapezoid
% 6 nu, corr and sigma

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
end
  
1; 

save('./fenics_inputs/inputs.mat','nu','mode_qoi','theta','delta_q')

% Extract number of degrees of freedom and nodal coordinates
frid = fopen(['./mesh/mesh_' num2str(coarse_level) '.txt'],'r');
mesh_info = fscanf(frid,'%d',[1,2]);
xy_coarse = fscanf(frid,'%g %g',[2 mesh_info(1)]); xy_coarse = xy_coarse';
if mode_qoi == 0
    coarse_grid = 6; 
elseif mode_qoi == 1
    coarse_grid = size(xy_coarse,1);
elseif mode_qoi == 2
    coarse_grid = size(xy_coarse,1);
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



%Uf fine
if mode_qoi == 0  
    load('L_data/Idx_f')
    load('L_data/Idx_c')
    load('L_data/Uf_line')
    Uf = U(Idx_f,:); 

%     load('L_data/Uf_line_2')
%     Uf = Uf(Idx_f,:); 
    
    Uc = Uc(Idx_c,:); 
elseif mode_qoi == 1
    load('L_data/Uf_field')
    Uf = U; 
elseif mode_qoi == 2
    load('L_data/Uf_stress')
    Uf = Uf; 
end


samps = 1:nsim; 

Uf = Uf(:,samps);
Uc = Uc(:,samps);

1; 
mean(Uc(:))


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

1; 

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

% U = Uf(:,ix_r)*P_s_r; 
% % save('Beam_design/Bi_nom','Bi')
% save('L_data/Ub_stress_opt','U')

err_Ahat
ahat_error_est/norm(A_R)
1; 



end