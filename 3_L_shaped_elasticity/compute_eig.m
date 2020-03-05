function [lambda,z,x] = compute_eig(sigma,cl,NTerms,NPoints)
%Returns eigenvalues of the correlation function and eigenfunctions at
%locations corresponding to the nodes of the Galerkin projection.

dx = 1/(NPoints-1);

C(1:NPoints,1:NPoints) = 0;
M(1:NPoints,1:NPoints) = 0;

x(1:NPoints) = dx*(0:(NPoints-1));

%Fill C
for i=1:NPoints
    i;
    for j=1:i
        C(i,j) = Cfunc(x(i),x(j),dx,sigma,cl);
        C(j,i) = C(i,j);
    end
end


%Fill M
for i=1:NPoints
    for j=1:i
        if (i==j)
            if (i==1 || i==NPoints)
                M(i,j) = 1/3*dx;
            else
                M(i,j) = 2/3*dx;
            end
        elseif (abs(i-j)==1)
            M(i,j) = 1/6*dx;
        end
        M(j,i) = M(i,j);
    end
end

%Get the eigenvalues and the eigenvectors of the Galerikin projection
[z,Lambda] = eigs(C,M,NTerms);

% Make the sign of all eigenvectors the same the same (needed for bi-fidelity construction)
for j=1:NTerms
    z(:,j) = z(:,j)/sign(z(1,j)); 
end


lambda = diag(Lambda);