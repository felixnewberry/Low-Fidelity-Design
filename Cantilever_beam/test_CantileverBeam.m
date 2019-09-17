close all;  clear all;  clc;

% load data matrices
% load ./CantileverBeam/Uf.mat; U_f = Uf;
% load ./CantileverBeam/Uc.mat; U_c = Uc;
load ./Beam_data/Uf.mat; U_f = Uf;
load ./Beam_data/Uc.mat; U_c = Uc;

%%%% changes made by Felix: 
% just first 100 samples
% matrixIDvR as opposed to matrixID
samps = 1:100; 
U_f = U_f(:,samps);
U_c = U_c(:,samps);I

% problem parameters
m = size(U_c,1);      % number of rows of A coarse
n = size(U_c,2);      % number of cols of A
r = min(m,n); % maximal rank of A
tol = 1e-3;   % approximate tolerance for matrix ID (it should give correct order of magnitude)
R = 1; 

% generate an approximately low-rank matrix
[U,S,V] = svd(U_c);
Sdiag = diag(S(1:r,1:r))/S(1,1);

% plot the singular values of A
semilogy(Sdiag,'ro')
title('Singular values of U_c')

% Show how matrix ID works
% [P,ix] = matrixID(U_c,tol^2);
[P,ix] = matrixIDvR(U_c,R);
Uc_id = U_c(:,ix) * P;
Uf_id = U_f(:,ix) * P;

err_id_c = norm(U_c - Uc_id,'fro')/norm(U_c,'fro');
err_id_f = norm(U_f - Uf_id,'fro')/norm(U_f,'fro');

err_id_f = norm(U_f - Uf_id)/norm(U_f)


fprintf('*** RESULTS ***\n')
fprintf(' Desired accuracy   = %e\n', tol)
fprintf(' Accuracy of ID for coarse model    = %e\n', err_id_c) 
fprintf(' Accuracy of ID for fine model    = %e\n', err_id_f) 
fprintf(' Approximation rank = %d\n', length(ix));


figure;
nc = size(U_c,1);
nf = size(U_f,1);

plot((0+1/(nc+1):1/(nc+1):1-1/(nc+1))',U_c(:,23:25),'b.-'); hold on;
plot((0+1/(nf+1):1/(nf+1):1-1/(nf+1))',U_f(:,23:25),'r-'); hold on;
plot((0+1/(nf+1):1/(nf+1):1-1/(nf+1))',Uf_id(:,23:25),'k--')


figure;
plot((0+1/(nc+1):1/(nc+1):1-1/(nc+1))',var(U_c'),'b.-'); hold on;
plot((0+1/(nf+1):1/(nf+1):1-1/(nf+1))',var(U_f'),'r')
plot((0+1/(nf+1):1/(nf+1):1-1/(nf+1))',var(Uf_id'),'k--')
