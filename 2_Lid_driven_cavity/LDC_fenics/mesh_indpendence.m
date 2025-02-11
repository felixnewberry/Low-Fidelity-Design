%%% Mesh indpependence

clear all
close all
clc
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
LW_axis = 1; 

% Colors
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; 
c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880]; 
c6 = [0.3010, 0.7450, 0.9330]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain data - edit file for different QoI or different mesh sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P Mid and P Base

pyFi_mesh = 'sudo python3.6 run_LDC_mesh_independence.py'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QoI_vec = 0:7;

% notes: u field - needs higher resolution - down to 1 % for 64 compare to
% 256 - do two more to check. 
% U field - 1 % at 64 comparing to 256 - could go higher... 
% P field - terrible - singularity in corners
% U mid - 0.09 % at 64
% U vert - x velocity? 
% P mid - 0.24 % at 64
% P vert - 0.2 % at 64... not necessarily a useful QoI... 
% P top - terrible - singularities 


% First - save high -fidelity for the two U, then for P vert and P base. 

% New set of QoI to consider: Save all of these. % set up. 
% U mid v
% U vert u
% P mid
% P base
% P vert

error_mat = zeros(8,6); 

for i_qoi = 1:length(QoI_vec)
    QoI = QoI_vec(i_qoi); 
    
if QoI == 0
    load('u_meshes/u_mesh_u_field')
    plot_label = 'U Field';
elseif QoI == 1
    load('u_meshes/u_mesh_P_field')
    plot_label = 'P Field';
elseif QoI == 2
    load('u_meshes/u_mesh_u_mid')
    plot_label = 'U Mid';
elseif QoI == 3
    load('u_meshes/u_mesh_P_mid')
    plot_label = 'P Mid';
elseif QoI == 4
    load('u_meshes/u_mesh_P_top')
    plot_label = 'P Top';
elseif QoI == 5
    load('u_meshes/u_mesh_P_vert')
    plot_label = 'P Vert';
elseif QoI == 6
    load('u_meshes/u_mesh_P_base')
    plot_label = 'P Base';
elseif QoI == 7
    load('u_meshes/u_mesh_u_vert')
    plot_label = 'u vert';
end

load('LDC_data/x_64')
x_64 = x_64(:,1); 

% plot u_matrix

nx_vec = [4,8,16,32,64,128,256]; 

% Calculate error
error = zeros(length(nx_vec)-1,1); 
for i = 1:length(nx_vec)-1
    error(i) = norm(u_matrix(i,:)-u_matrix(end,:))/norm(u_matrix(end,:));
end

figure
semilogy(nx_vec(1:end-1),100*error,'-s','Color',c2,'LineWidth',LW);
ylabel('Relative Error $[\%]$', 'interpreter', 'latex', 'fontsize', FS)
xlabel('$nx$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on; 
axis tight;
set(gcf,'Position',size_1)
title(strcat(plot_label,' Convergence'),'Interpreter','latex')
    
results = [nx_vec(1:end-1);100*error']'

error_mat(i_qoi,:) = 100*error'; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate high-fidelity data: first save x_64... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% larger distribution of random variables 

% prevous nu_vec and u_top_vec were 5 % and 10 % 
% try 10 % and 20 % to mix things up. 

n_samples = 200; 

nu_nom = 0.01; 
u_nom = 1; 

delta_nu = 0.1; 
delta_u = 0.2; 

nu_vec = nu_nom*(1+delta_nu*(rand(n_samples,1)*2-1)); 
u_top_vec = u_nom*(1+delta_u*(rand(n_samples,1)*2-1)); 

% save('LDC_data/u_nu_vec_2','nu_vec', 'u_top_vec')

figure
p1 = semilogy(nx_vec(1:end-1),error_mat(1,:),'-s','Color',c1,'LineWidth',LW);
hold on;
p2 = semilogy(nx_vec(1:end-1),error_mat(2,:),'--s','Color',c2,'LineWidth',LW);
p3 = semilogy(nx_vec(1:end-1),error_mat(3,:),':o','Color',c3,'LineWidth',LW);
p4 = semilogy(nx_vec(1:end-1),error_mat(4,:),'-.d','Color',c4,'LineWidth',LW);
p5 = semilogy(nx_vec(1:end-1),error_mat(5,:),':v','Color',c5,'LineWidth',LW);
p6 = semilogy(nx_vec(1:end-1),error_mat(6,:),'--<','Color',c6,'LineWidth',LW);
p7 = semilogy(nx_vec(1:end-1),error_mat(7,:),'-.>','LineWidth',LW);
p8 = semilogy(nx_vec(1:end-1),error_mat(8,:),':*','LineWidth',LW);
hold off; 
ylabel('Relative Error $[\%]$', 'interpreter', 'latex', 'fontsize', FS)
xlabel('$nx$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on; 
axis tight;
set(gcf,'Position',size_1)
legend([p1,p2,p3,p4,p5,p6,p7,p8],{'U Field','P Field','U Mid','P Mid','P Top', 'P Vert', 'P Base', 'U Vert'},...
       'interpreter', 'latex', 'fontsize', FS_leg/2)
% title(strcat(plot_label,' Convergence'),'Interpreter','latex')

