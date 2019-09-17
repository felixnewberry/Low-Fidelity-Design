function [error_bound,err_Ahat,efficacy] = my_beam_bound(X,nsim, n, r, mode,n_reps)

%%% Inputs

% X - vector of deltas 
% X(1) = delta t1= delta t2; 
% X(2) = delta t3; 
% nsim - number of runs within one sample, ie 100 for cantilever beam
% n - truncation - subset of data from which bound is estimated. 
% r - truncation with which bi-fidelity model is found

% mode:
% 0 test w
% 1 test t1
% 2 test t2 
% 3 test t3
% 4 test 

% n_reps number of repetitions of bound samples

%%% Outputs
% error_bound - average over repetitions
% err_Ahat - average over repetitions
% efficiacy - average over repetionns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beam details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unifom load 
q0 = 10;
s4 = 1;

% Generate beam samples with new inputs
% Original sizes

w = 1; % Cross section width
t1 = 0.2; % height of part 1
t2 = 0.2;
t3 = 5;
    
if mode == 0
    w = w*(1+X);
elseif mode == 1
    t1 = t1*(1+X);
elseif mode == 2
    t2 = t2*(1+X);
elseif mode == 3
    t3 = t3*(1+X);
elseif mode == 4
    t1 = t1*(1+X(1));
    t2 = t2*(1+X(1));
    t3 = t3*(1+X(2));
end

L = 50;

% Young modulus E_j = E0_j + s_j * Z_j and Z_j~U[-1,1]
E01 = 1e6; 
s1 = 1e5;
E02 = 1e6; 
s2 = 1e5;
E03 = 1e4; 
s3 = 1e3;

% coordinates
load('Beam_data/x_highfidelity.txt')

%Uf fine
load('Beam_data/Uf')

% xi 
load('Beam_data/xi')

% xi = load('/Users/felixnewberry/Google Drive/1_Research/3_low_fidelity_design/CantileverBeam/Beam_data/xi.mat');
% xi = xi.xi; 
% 
% load('/Users/felixnewberry/Google Drive/1_Research/3_low_fidelity_design/CantileverBeam/Beam_data/Uf.mat','Uf')
% load('/Users/felixnewberry/Google Drive/1_Research/3_low_fidelity_design/CantileverBeam/x_highfidelity.txt','x_highfidelity')


for isim = 1:nsim
    E1 = E01 + s1 * xi(isim,1);
    E2 = E02 + s2 * xi(isim,2);
    E3 = E03 + s3 * xi(isim,3);
    q  = q0  + s4 * xi(isim,4);
    Uc(:,isim) = EB_Cantilever(L,t1,t2,t3,w,E1,E2,E3,q,x_highfidelity);
end

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