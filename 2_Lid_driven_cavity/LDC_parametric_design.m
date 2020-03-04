clear all
close all
clc

% Objective:
% Minimize epsilon tau over low-fidelity parameterization of model. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choose optimization method and mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% is this way out of date?? 

% record Bound:bi:low
% 0 is u mid v 
% 1 is u vert u 
% 2 is P mid
% 3 is P base           
% 4 is P vert

% convert point_test into nom_opt
individual_search = 0; % 1 for nu, 2 for u
grid_search = 0; % 
nom_opt = 1; % Save data for nominal and optimal runs - have to edit bound too

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;     % Line width
MS = 8;     % Marker Size
FS_leg = 16; % Font size legend


size_1 = [0,0,670,515]; 
size_2 = [0,0,1340,515]; 

FS = 28;    % Font size axis
FS_axis = 18; 
LW_axis = 2; 

% Colors
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; 
c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880]; 
c6 = [0.3010, 0.7450, 0.9330]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate samples 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 200 samples

% % n_samples = 200; 
% % % rho = 997; 
% % % mu_av = rho/Re; 

% Re = 100; 
% nu_av = 1/Re; 
% u_top_av = 1; 

% xi_mat = (rand(n_samples,2)-0.5)*2; 
% nu_vec = nu_av*(1+0.05*xi(:,1)); 
% u_top_vec = u_top_av*(1+0.1*xi(:,2)); 

% save('nu_vec','nu_vec')
% save('u_top_vec','u_top_vec')
% save('xi_mat','xi_mat')

% % nx = 32 took 133.21 s to run. 
% % nx = 8 took 7.132 s to run. 
% % nx = 6 took 4.94 s to run. 
% % nx = 4 took 3.88 s to run. 

% % % load('LDC_data/nu_vec.mat')
% % % load('LDC_data/u_top_vec.mat')
% % % load('LDC_data/xi_mat.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data, vanilla approach first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordinates x_32 and y_32
load('LDC_data/x_32.mat');
load('LDC_data/y_32.mat');
load('LDC_data/x_all.mat');

% random inputs, xi_mat, nu_vec, u_top_vec
% load('LDC_data/xi_mat.mat');
load('LDC_data/u_nu_vec.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Individual Sensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if individual_search >= 1
    
    % Available hyper-parameters are varying nu and u_top. 
    % start with variation +- 10 
    
    nx = 4; 
%     r = 1; 
    r = 1; 
    n = r+2; 
        
    
%     delta_u_vec = -0.6:0.02:0.0;
%     delta_nu_vec = 0.0:0.1:3;
    
    % broke: %     p1 = semilogy(1:length(s),s,'-ob','LineWidth',LW);
    
    if individual_search == 1
        delta_nu_vec = -0.7:0.2:3;
%         delta_nu_vec_1 = delta_nu_vec_0; 
        delta_u_vec = zeros(length(delta_nu_vec), 1); 
    elseif individual_search == 2
       delta_u_vec = -0.8:0.2:2.4;  
        delta_nu_vec = zeros(length(delta_u_vec), 1); 
%         delta_nu_vec_1 = delta_nu_vec_0; 
    elseif individual_search == 3
% %         delta_nu_vec_0 = -0.5:0.1:0.5;
%         delta_nu_vec_0 = -0.9:0.1:3;
% 
%         delta_nu_vec_1 = zeros(length(delta_nu_vec_0), 1); 
%         delta_u_vec = delta_nu_vec_1; 
    elseif individual_search == 4
% %         delta_nu_vec_0 = -0.5:0.1:0.5;
%         delta_nu_vec_1 = -0.9:0.1:3;
% 
%         delta_nu_vec_0 = zeros(length(delta_nu_vec_1), 1); 
%         delta_u_vec = delta_nu_vec_1; 
    end
    
    n_samps = length(delta_nu_vec); 
    
%     save('LDC_design/delta_nu_vec','delta_nu_vec')
%     save('LDC_design/delta_u_vec','delta_u_vec')
%     save('LDC_design/delta_P_mid_nu_vec','delta_nu_vec')
%     save('LDC_design/delta_P_mid_u_vec','delta_u_vec')
%     
    error_bound_mat = zeros(5,n_samps); 
    error_Bi_mat = zeros(5,n_samps);    
    
%     run_count = 1; 
    
    for i_param = 1:n_samps

        [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r, delta_u_vec(i_param), delta_nu_vec(i_param), nom_opt);
        
%         run_count = run_count+1; 
        
        error_bound_mat(:,i_param) = error_bound; 
        error_Bi_mat(:,i_param) = err_bi; 
        
        fprintf("Percent complete: %d \n",i_param/n_samps*100);       
    end
%     
 1; 
 
if individual_search == 1
    plot_label = '$ \Delta \nu [\%]$';
    delta_vec = delta_nu_vec; 
%     save('LDC_design/line_qoi_nu_1','error_bound_mat', 'delta_vec')
    save('LDC_design/individual_qoi_nu','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_nu_2_n8','error_bound_mat', 'delta_vec')

elseif individual_search == 2
    plot_label = '$ \Delta U [\%]$';
    delta_vec = delta_u_vec; 
%     save('LDC_design/line_qoi_u_1','error_bound_mat', 'delta_vec')
    save('LDC_design/individual_qoi_u','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_u_2_n8','error_bound_mat', 'delta_vec')
elseif individual_search == 3
%     plot_label = '$ \Delta \nu_0 [\%]$';
%     delta_vec = delta_nu_vec_0; 
% %     save('LDC_design/line_qoi_u_1','error_bound_mat', 'delta_vec')
% %     save('LDC_design/line_qoi_u_2','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_nu_y_n4_2','error_bound_mat', 'delta_vec')
elseif individual_search == 4
%     plot_label = '$ \Delta \nu_1 [\%]$';
%     delta_vec = delta_nu_vec_1; 
% %     save('LDC_design/line_qoi_u_1','error_bound_mat', 'delta_vec')
% %     save('LDC_design/line_qoi_u_2','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_nu_y1_n4_2','error_bound_mat', 'delta_vec')
end

1; 
% Plot error bound 
figure
hold on
p1 = plot(100*delta_vec,100*error_bound_mat(1,:), 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(100*delta_vec,100*error_bound_mat(2,:),'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_bound_mat(3,:), 'Color',c3,'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_vec,100*error_bound_mat(4,:), 'Color',c4,'LineWidth',LW,'MarkerSize',MS); 
p5 = plot(100*delta_vec,100*error_bound_mat(5,:), 'Color',c5,'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Base', '$P$ Vert'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
grid on


% Plot bi-fidelity error 
figure
hold on
p1 = plot(100*delta_vec,100*error_Bi_mat(1,:), 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(100*delta_vec,100*error_Bi_mat(2,:),'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_Bi_mat(3,:), 'Color',c3,'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_vec,100*error_Bi_mat(4,:), 'Color',c4,'LineWidth',LW,'MarkerSize',MS); 
p5 = plot(100*delta_vec,100*error_Bi_mat(5,:), 'Color',c5,'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bi $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Base', '$P$ Vert'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
grid on


%     s = svd(Uc); s = s/s(1); 
%     
%   % Plot svd: 
%     figure
%     p1 = semilogy(1:length(s),s,'-ob','LineWidth',LW);
%     xlabel('index $i$', 'interpreter', 'latex', 'fontsize', FS)
%     ylabel('Normalized Singular Value', 'interpreter', 'latex', 'fontsize', FS)
%     grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
%     
    % 200 samples r = 10. 
    % err_Ahat nx = 8 = 1.55% % error_low is 2.72 % 1.20821 s
    % err_Ahat nx = 6 = 1.58% % error_low is 5.94 
    % err_Ahat nx = 4 = 1.67% % error_low is 18.58 
    
    % LINE PLOT TAKE AWAYS... 
end

if grid_search == 1
    
    nx = 4; 
    r = 1; 
    n = r+2; 
    
    u_lim = [-0.8,2.4]; 
    nu_lim = [-0.7,3.0];
    
    n_grid = 20; 
    
    delta_nu_vec = linspace(nu_lim(1),nu_lim(2),n_grid); 
    delta_u_vec = linspace(u_lim(1),u_lim(2),n_grid); 
    
%     delta_u_vec = zeros(1, length(delta_nu_vec)); 
%     delta_nu_vec_0 = delta_nu_vec; 
%     delta_nu_vec_1 = delta_nu_vec; 
%     
    error_bound_mat = zeros(5,length(delta_nu_vec),length(delta_u_vec)); 
    error_Bi_mat = zeros(5,length(delta_nu_vec),length(delta_u_vec)); 
    
%     for i_nu = 1:length(delta_nu_vec)
%         for i_u = 1:length(delta_u_vec)
%     run_count = 1; 
    
    for i_nu = 1:length(delta_nu_vec)
        for i_u = 1:length(delta_u_vec)

%             delta_u_vec(i_u)
%             delta_nu_vec(i_nu)
    
            [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r,delta_u_vec(i_u),delta_nu_vec(i_nu),nom_opt);
            error_bound_mat(:,i_nu, i_u) = error_bound; 
            error_Bi_mat(:,i_nu, i_u) = err_bi; 
            
%             run_count = run_count+1; 
            1; 
        
            
        end 
        fprintf("Percent complete: %d \n",i_nu/length(delta_nu_vec)*100);  
    end
    
%         save('LDC_design/grid_search_2_n8','error_bound_mat', 'error_Bi_mat', 'delta_u_vec','delta_nu_vec')
%         save('LDC_design/grid_search_test','error_bound_mat', 'error_Bi_mat', 'delta_u_vec','delta_nu_vec')
%         save('LDC_design/grid_search_nu_linear_test','error_bound_mat', 'error_Bi_mat', 'delta_u_vec','delta_nu_vec_0','delta_nu_vec_1')
%         save('LDC_design/grid_search_2_Transform','error_bound_mat', 'error_Bi_mat', 'delta_u_vec','delta_nu_vec')
        save('LDC_design/grid_search','error_bound_mat', 'error_Bi_mat', 'delta_u_vec','delta_nu_vec')

end

if nom_opt == 1
    
nx = 4; 
r = 1; 
n = r+2; 

% Find nominal and optimal data - just for U mid and P base 
delta_u_vec = [0, 1.5579, 0.71579]; 
delta_nu_vec = [0, -0.11579, 2.8053]; 

error_bound_mat = zeros(length(delta_u_vec),1); 
error_Bi_mat = zeros(length(delta_u_vec),1); 

nom_opt_count = nom_opt; 

for i_test = 1:length(delta_u_vec)
     
[error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r,delta_u_vec(i_test),delta_nu_vec(i_test),nom_opt_count);

nom_opt_count = nom_opt_count+1; 
% 
% error_bound_mat(i_test) = error_bound; 
% error_Bi_mat(i_test) =err_bi; 
end

% error_bound_mat
% error_Bi_mat
end




