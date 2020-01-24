clear all 
close all
clc

%%% Test the error bound, why is it failing?

rand_sample = 1:n; 

% % Normalize matrices
% Uc= Uc;
% Uf = Uf;

% Bi-fid matrices. 
B_R = Uc(:,rand_sample)/norm(Uc,'fro');
A_R = Uf(:,rand_sample)/norm(Uf,'fro');