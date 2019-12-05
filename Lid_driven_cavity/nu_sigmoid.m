clear all
close all
clc

%  nu test
nu_orig = 0.01; 
delta_nu_0 = 0; 
delta_nu_1 = -0.5; 

nu_0 = nu_orig*(1+delta_nu_0); 
nu_1 = nu_orig*(1+delta_nu_1); 

x_vec = linspace(0,1,100);
[X,Y] = meshgrid(x_vec); 

nu_mat_1 = zeros(length(x_vec), length(x_vec)); 

for i_x = 1:length(x_vec)
    for i_y = 1:length(x_vec)
        nu_mat_1(i_x,i_y) = nu_0+(nu_1-nu_0)*x_vec(i_y);
    end
end

figure
contourf(X,Y,nu_mat_1')
colorbar

s_1 = -1./(1+exp(-(x_vec)));
figure
plot(x_vec, s_1)

nu_mat_2 = zeros(length(x_vec), length(x_vec)); 

s_tight = 10; 
s_center = 0.5; 
vort = [0.6215,0.7357]; 

for i_x = 1:length(x_vec)
    for i_y = 1:length(x_vec)
        d_1 = sqrt((x_vec(i_x) - vort(1))^2+(x_vec(i_y) - vort(2))^2);
        nu_mat_2(i_x,i_y) = nu_0*(1+delta_nu_1/(1+exp(s_tight*(d_1-s_center))));
    end
end

figure
contourf(X,Y,nu_mat_2')
colorbar


x_vec = 0:0.01:1;
s_tight = 10; 
s_center = 0.5; 

s_1 = 1./(1+exp(s_tight*(x_vec-s_center)));
figure
plot(x_vec, s_1)

x_vec = 0:0.01:1;
s_tight = 10; 
s_center = 0.5; 

s_1 = 1./(1+exp(s_tight*(x_vec-s_center)));
figure
plot(x_vec, s_1)

% next try a ring... 