function youngs_samp_gen_gaussian_cov(youngs_0, sigma, L_c, d, xy_coord,...
    xi,num_kl_grid,delta)

% note code applied for square domain but will read youngs over the l-shaped

%   Variables:
%   L_c = correlation length of the process
%   sigma = sqrt of standard deviation of the temperature
%   xi = a N by d^2 matrix of N independent realization of d KL random
%   variables
%   num_kl_grid = Number of grid point for computing kL examsion (not related to the PDE
% grid)

% Outputs:
%   youngs_samp = a matrix containing the samples of youngs

% Number of terms in expansion
NTerms = d*d;

% Number of samples
[nsim, ~] = size(xi);

% Generate the grid along y direction
x = xy_coord(:,1);
y = xy_coord(:,2);

% Number of grid points in each direction
nx = length(x);

% First generate the KL eigen-modes of D=1

[lambda_1d,evec_1d,kl_grid] = compute_eig(sqrt(sigma),L_c,d,num_kl_grid);

% Initialize matrices and interpolate the eigenvectors
evec_x = zeros(nx,d);
evec_y = zeros(nx,d);

for i=1:d
    evec_x(:,i) = interp1(kl_grid,evec_1d(:,i),(x+1)/2);
    evec_y(:,i) = interp1(kl_grid,evec_1d(:,i),(y+1)/2);
end


% Eigenvalues for x and y d
lambda_xy = zeros(d,d); % Generates the tensor product of eigenvalues
for i=1:d
    for j=1:d
        lambda_xy(i,j) = lambda_1d(i,1)*lambda_1d(j,1);
    end
end

% Rearrange \lmbda_i*\lmbda_j in a descend order (not needed when NTerm =
% d^2)
lambda = zeros(NTerms,3);
[lambda(:,1), idx] = sort(lambda_xy(:),'descend');
for i=1:NTerms
    [I, J] = ind2sub(size(lambda_xy),idx(i));
    lambda(i,2:3) = [I,J];  %gives the index of the evec
end

% Generate the deterministic part of the KL expansion of the gaussian field
g = zeros(nx,NTerms);
for j=1:NTerms
    for i=1:nx
        g(i,j) = sqrt(lambda(j,1))*evec_y(i,lambda(j,3))*evec_x(i,lambda(j,2));
    end
end

%save data.mat lambda lambda_1d evec_1d g
      
% Generate the realizations of the lognormal youngs for FENICS
Youngs = zeros(nsim,nx);

for i=1:nsim
    youngs = zeros(1,nx);
    % Next loop generates the realizations of the Gaussian
    for k=1:NTerms
        youngs = youngs + g(:,k)'*xi(i,k);
    end
    % Genrate the lognormal youngs
    Youngs(i,:) = youngs_0 + exp(youngs);
end

Youngs = Youngs*(1+delta); 
1; 
save('./fenics_inputs/Youngs.mat', 'Youngs')

end



