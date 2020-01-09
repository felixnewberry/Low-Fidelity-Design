clear all
close all
clc

% Objective: Process low-fidelity data to explore the transformation matrix
% projection.
% I should also test the transformation vs interpolative decomposition. -
% it outperforms decomposition - but lacks the bound. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load and process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% U_matrix 0 through 7 corresponds to u_y, u_x, p_mid, p_vert, p_base,
% p_field, u_x_field, u_y_field

% Step through samples
n_samples = 200; 
n_run = 400; 

n_x = 65; 
n_xy = 4225; 

u_y = zeros(n_samples, n_x, n_run); 
u_x = zeros(n_samples, n_x, n_run); 
p_mid = zeros(n_samples, n_x, n_run); 
p_vert = zeros(n_samples, n_x, n_run); 
p_base = zeros(n_samples, n_x, n_run); 
p_field = zeros(n_samples, n_xy, n_run); 
u_x_field = zeros(n_samples, n_xy, n_run); 
u_y_field = zeros(n_samples, n_xy, n_run); 


for i_sample = 1:n_run
    i_sample
    load(strcat('u_meshes/L_data/u_matrix_[[',num2str(1),']].mat'))
    u_y(:,:,i_sample) = u_matrix_0; 
    u_x(:,:,i_sample) = u_matrix_1; 
    p_mid(:,:,i_sample) = u_matrix_2; 
    p_vert(:,:,i_sample) = u_matrix_3; 
    p_base(:,:,i_sample) = u_matrix_4; 
    p_field(:,:,i_sample) = u_matrix_5; 
    u_x_field(:,:,i_sample) = u_matrix_6; 
    u_y_field(:,:,i_sample) = u_matrix_7; 
end


%%% Save data
save('u_meshes/L_mat/u_y', 'u_y')
save('u_meshes/L_mat/u_x', 'u_x')
save('u_meshes/L_mat/p_mid', 'p_mid')
save('u_meshes/L_mat/p_vert', 'p_vert')
save('u_meshes/L_mat/p_base', 'p_base')

save('u_meshes/L_mat/p_field', 'p_field', '-v7.3')
save('u_meshes/L_mat/u_x_field', 'u_x_field', '-v7.3')
save('u_meshes/L_mat/u_y_field', 'u_y_field', '-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u_lim = [-0.8,2.4]; 
    nu_lim = [-0.7,3.0];
    
    n_grid = 20; 
    
    delta_nu_vec = linspace(nu_lim(1),nu_lim(2),n_grid); 
    delta_u_vec = linspace(u_lim(1),u_lim(2),n_grid); 
%     for i_nu = 1:length(delta_nu_vec_0)
%         for i_u = 1:length(delta_u_vec)