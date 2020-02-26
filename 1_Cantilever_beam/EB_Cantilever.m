function [disp] = EB_Cantilever(L,t1,t2,t3,w,E1,E2,E3,q,x)

% Calculate the young's modulus ratios and equivalen areas (now all section
% are E3)

n1 = E1/E3;
n2 = E2/E3;

w1 = w * n1;
w2 = w * n2;

% Calculate the position of neutral axis

area = t1 * w1 + t2 * w2 + t3 * w; 
yb = (t1 * w1)*(t1/2 + t2 + t3) + (t2 * w2)*(t2/2) + (t3 * w)*(t3/2 + t2);
yb = yb/area;

% Compute the moment of inertia I
I = (1/12 * w1 * t1^3)+ (t1 * w1)*(t1/2 + t2 + t3 - yb)^2 + ...
    (1/12 * w2 * t2^3)+ (t2 * w2)*(t2/2 - yb)^2 + ...
    (1/12 * w  * t3^3)+ (t3 * w )*(t3/2 + t2 - yb)^2;

% could parameterize I. 

disp = (-q*x.^4/24 + q*L*x.^3/6 - q*L^2*x.^2/4)/(E3*I);



end

