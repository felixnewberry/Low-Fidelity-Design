function [Tau, S_1] = get_sobol_indices(c, index_pc)

% Inputs:
% c - PCE coefficient vector
% index_pc -  pce basis choice

% Outputs:
% Tau - total index
% S_1 - First order index

%c = Psi\U;

var_f = 0;


D = size(index_pc,2); % stochastic dimension
P = size(index_pc,1); % number of terms in pc expansion

% Calculate variance
for j = 2:P
    var_f = var_f+  c(j)^2; %*mean(Psi(:,j).^2); %var_f + (U(j)- mu_f)^2;
end
%var_f = var_f/D;

Tau = zeros(D,1);
S_1 = Tau;

for d = 1:D
   coeff_inds = find(index_pc(:,d) ~= 0);
   for j = 1:size(coeff_inds)
        Tau(d) = Tau(d) + c(coeff_inds(j))^2;
   end
   
end

for d = 1:D
    %coeff_inds = []
    for n = 1:P 
       if index_pc(n,d) ~= 0 && sum(index_pc(n,setdiff(1:D,d))) == 0
            S_1(d) = S_1(d) + c(n)^2;        
            %coeff_inds = [coeff_inds n];
       end
    end
    
    %coeff_inds
end

Tau = Tau./var_f
S_1 = S_1./var_f

end

