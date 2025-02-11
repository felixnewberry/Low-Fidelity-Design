function [error_bound,err_Ahat,efficacy] = my_beam_bound(X,nsim, n, r, mode)

%%% Inputs

% X - vector of deltas 
% X(1) = delta h1= delta h2; 
% X(2) = delta h3; 
% nsim - number of runs within one sample, ie 100 for cantilever beam
% n - truncation - subset of data from which bound is estimated. 
% r - truncation with which bi-fidelity model is found

% mode:
% 0 test w
% 1 test h1
% 2 test h2 
% 3 test h3
% 4 test h1 = h2
% 5 test h1 = h2 and h3
% 6 save nominal and optimal results

%%% Outputs
% error_bound 
% err_Ahat 
% efficiacy 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beam details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unifom load 
q0 = 10;
s4 = 1;

% Generate beam samples with new inputs
% Original sizes

w = 1; % Cross section width
h1 = 0.2; % height of part 1
h2 = 0.2;
h3 = 5;
    
if mode == 0
    w = w*(1+X);
elseif mode == 1
    h1 = h1*(1+X);
elseif mode == 2
    h2 = h2*(1+X);
elseif mode == 3
    h3 = h3*(1+X);
elseif mode == 4
    h1 = h1*(1+X(1));
    h2 = h2*(1+X(1));
elseif mode == 5
    h1 = h1*(1+X(1));
    h2 = h2*(1+X(1));
    h3 = h3*(1+X(2));
elseif mode == 6        % special case for running nom and opt - save data at end of script
    h1 = h1*(1+X(1));
    h2 = h2*(1+X(1));
    h3 = h3*(1+X(2));
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
load('Beam_data/x_highfidelity.txt','x_highfidelity')

%Uf fine
load('Beam_data/Uf','Uf')

% xi 
load('Beam_data/xi','xi')

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
    Uc(:,isim) = EB_Cantilever(L,h1,h2,h3,w,E1,E2,E3,q,x_highfidelity);
end

samps = 1:nsim; 

Uf = Uf(:,samps);
Uc = Uc(:,samps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fid details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

% Obtain column skeleton of P
[P_s,ix] = matrixIDvR(B,r);

% Error bound inputs
normC = norm(P_s);
sb = svd(B); 
err_Bhat = norm(B-B(:,ix)*P_s); 
N = nsim;

% Subset of vectors for bi-fidelity error estimate 
% Ensure column skeleton slection are used + additional samples to reach
% total of n
% rand_sample = [ix, getfield(setxor(ix,1:N), {1:n-numel(ix)})];
rand_sample = [ix, getfield(setxor(ix,1:n), {1:n-numel(ix)})];

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

% % % Compute epsilon tau... 
[~, ahat_error_est,~, ~,~] = ...
    mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);

error_bound = ahat_error_est/norm(A);

% Do bi-fidelity estimate to compare with bound
err_Ahat = norm(A-A(:,ix)*P_s)/norm(A);
efficacy = error_bound/err_Ahat;

1; 

if mode == 6
    Ub = Uf(:,ix)*P_s;
  
    % relative error of tip deflection
    errors_Ahat_1 = (A-A(:,ix)*P_s);
    tip_error_bi = errors_Ahat_1(end,:)./A(end,:);
    
    errors_L = (Uf-Uc);
    tip_error_L = errors_L(end,:)./Uf(end,:);
    
    if X(1) == 0 
        save('Beam_design/Nom','Uc', 'Ub', 'sb','rand_sample','n','r')
       	save('Beam_design/tip_error_nom','tip_error_bi','tip_error_L','rand_sample','n','r')
    else
        save('Beam_design/Opt','Uc', 'Ub', 'sb','X','rand_sample','n','r')
        save('Beam_design/tip_error_opt','tip_error_bi','tip_error_L','X','rand_sample','n','r')
    end
end

end
