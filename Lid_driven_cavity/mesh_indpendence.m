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

load('u_meshes/u_mesh3')
load('LDC_data/x_32')
x_32 = x_32(:,1); 

% plot u_matrix

nx_vec = [4,8,16,32,64,128,256]; 

figure
p1 = plot(x_32,u_matrix(1,:),'Color',c1,'LineWidth',LW);
hold on
p2 = plot(x_32,u_matrix(2,:),'Color',c2,'LineWidth',LW);
p3 = plot(x_32,u_matrix(3,:),'Color',c3,'LineWidth',LW);
p4 = plot(x_32,u_matrix(4,:),'Color',c4,'LineWidth',LW);
p5 = plot(x_32,u_matrix(5,:),'Color',c5,'LineWidth',LW);
p6 = plot(x_32,u_matrix(6,:),'Color',c6,'LineWidth',LW);
p7 = plot(x_32,u_matrix(7,:),'LineWidth',LW);

ylabel('$P$ Base', 'interpreter', 'latex', 'fontsize', FS)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on; 
axis tight;


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

    
results = [nx_vec(1:end-1);100*error']'; 

% use 64 for high. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate high-fidelity data: first save x_64... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
