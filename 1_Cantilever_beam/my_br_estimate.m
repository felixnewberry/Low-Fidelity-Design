function [error_low,error_bi, u_bi] = my_br_estimate(r, N_hi, sigma, Uc, Uf, psi_ref, c_low, ix)

% Which samples to use? 
% ix + remainder of n used in bound + remainder ... 

% sample = datasample(1:size(Uf,1), N_hi, 'Replace', false); 
sample = [ix, getfield(setxor(ix,1:size(Uf,1)), {1:N_hi-numel(ix)})];

%             high fidelity model evaluations at sample
u_hi = Uf(sample,:); 

%             measurement matrix at sample - used to form high fid PCE
psi_hi= psi_ref(sample,:); 
 
% generate bi-fidelty estimate based on based on many low-fidelity samples
%             and some high fidelity. This is the general approach
[c_bi,alpha]= BR_FN(u_hi,psi_hi,c_low,r,sigma);

% This BR_FN function really really needs to be cleand up. 

% For given r and N combo, calculate errors
1;

psi_bi = psi_ref(:,2:end)*alpha';

% New reduced basis including column of ones
psi_bi = [ones(size(Uc, 1), 1) psi_bi]; 
% remove r+1 column before solving for coefficients.
psi_bi = psi_bi(:,1:end-1); 

% Multiply to find u_bi
u_bi = c_bi*psi_bi'; 
u_bi = u_bi'; 

1; 

%error_low = norm(Uc - Uf)/norm(Uf); 
error_low = 1; 
error_bi = norm(u_bi - Uf)/norm(Uf); 
end

