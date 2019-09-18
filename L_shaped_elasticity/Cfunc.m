function value = Cfunc(x1_0,x2_0,dx,sigma,cl)
%Returns the value of an element of the C matrix.

%Area where the linear finite elements are nonzero; everywhere else the
%integral is zero.
down1 = max(0,x1_0-dx);
up1 = min(1,x1_0+dx);

down2 = max(0,x2_0-dx);
up2 = min(1,x2_0+dx);

value = dblquad(@(x1,x2)sigma^2*exp(-(x1-x2).^2/cl.^2).*...
            (1-abs(x1-x1_0)/dx).*(1-abs(x2-x2_0)/dx),...
            down1,up1,down2,up2,1e-10);