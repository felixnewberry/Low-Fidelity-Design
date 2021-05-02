clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Cp for each ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eventually - Loop over ensembles

% folder_name_loop = 'limit_tests_pm5'; 
% folder_name_loop = 'limit_tests_pm10'; 
% folder_name_loop = 'batch_15_20'; 
% folder_name_loop = 'batch_50'; 
% folder_name_loop = 'batch_1000_all'; 
folder_name_loop = 'batch_500_airfoil'; 



% density, m, p, t, a, sa1 - 7
design_params_loop = importdata(strcat(folder_name_loop, '/design_params.dat'), ',' , 1); 
design_params_loop = design_params_loop.data; 

params_500_loop = importdata(strcat('limit_tests_pm5', '/params_500_orig.dat'), ',' , 1); 
params_500_loop = params_500_loop.data; 


% Set up for results

n_sims = size(params_500_loop,1); 
n_ens_loop = size(design_params_loop, 1); 
n_points = 200; 
xx = linspace(-1,1,n_points); 

cp_results = zeros(n_points, n_sims, n_ens_loop); 


for i_ens = 1:n_ens_loop
    
    % For ensemble: 
    load(strcat(folder_name_loop, '/airfoil_p_data_', num2str(i_ens),'.mat'))
    load(strcat(folder_name_loop, '/airfoil_x_data_', num2str(i_ens),'.mat'))
    load(strcat(folder_name_loop, '/airfoil_y_data_', num2str(i_ens),'.mat'))
    load(strcat(folder_name_loop, '/scalar_data_', num2str(i_ens),'.mat'))

    % adjust naca_parms for ensemble
    naca_params = params_500_loop.*(1+[design_params_loop(i_ens,2:5), 0]); 

    % Sort p data, compute cp and interpolate into 200 spaced points -1 to 1
    [ cp_results(:,:,i_ens)] = ...
        sort_airfoil( naca_params, airfoil_p_data, airfoil_x_data, airfoil_y_data, scalar_data); 

end

% % save('results_data/cp_results_lims_pm5', 'cp_results')
% % save('results_data/cp_results_lims_pm10', 'cp_results')
% % save('results_data/cp_results_lims_pm15_20', 'cp_results')
% % save('results_data/cp_results_lims_pm25_30', 'cp_results')
save('results_data/cp_results_500', 'cp_results')