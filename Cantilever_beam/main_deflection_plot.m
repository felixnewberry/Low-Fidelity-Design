clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LW = 2;     % Line width
FS_leg = 16; % Font size legend


size_1 = [0,0,670,515]; 
size_2 = [0,0,1340,515]; 

% size_1 = [0, 0, 500, 350]; 

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
%%% Generate a new Uc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Compute the deflection of a cantilever beam with composit cross section,
% % where the Young's modulus of the each section is random
% 
% % Sizes of beam cross section: 1 for top flang, 2 for bottom flang and 3
% % for the web
% 
% % % Original sizes
% % t1 = 0.2;
% % t2 = 0.2;
% % t3 = 5;
% % w = 1; % Cross section width
% 
% % adjusted
% t1 = 0.2150;
% t2 = 0.2150;
% t3 = 75;
% w = 1; % Cross section width
% 
% 
% % parameters it would be good to vary: t3 (down from 5, ie to 4)
% % w from 1 to 0.5 
% 
% % Young modulus E_j = E0_j + s_j * Z_j and Z_j~U[-1,1]
% E01 = 1e6; 
% s1 = 1e5;
% E02 = 1e6; 
% s2 = 1e5;
% E03 = 1e4; 
% s3 = 1e3;
% 
% % Length of the beam 
% L = 50;
% load('Beam_data/x_highfidelity.txt')
% load('Beam_data/xi')
% 
% %x = linspace(0,L,100)';
% 
% % Unifom load 
% q0 = 10;
% s4 = 1;
% 
% % Number of samples
% nsim = 100;
% 
% % Start the simulation
% % rand('state',3);
% % xi = 2*rand(nsim,4) - 1;
% 
% 
% for isim = 1:nsim
%     E1 = E01 + s1 * xi(isim,1);
%     E2 = E02 + s2 * xi(isim,2);
%     E3 = E03 + s3 * xi(isim,3);
%     q  = q0  + s4 * xi(isim,4);
%     Uc(:,isim) = EB_Cantilever(L,t1,t2,t3,w,E1,E2,E3,q,x_highfidelity);
% end
% 
% save('Beam_design/Uc_opt','Uc');

% Uc_plot = Uc; 
% 


1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Beam_data/x_highfidelity.txt')
load('Beam_data/xi')

load('Beam_design/Uc_opt')
Uc_trial = Uc;

%Uf? fine
load('Beam_data/Uf')
Uf = Uf(:,1:100);
% Uc course
load('Beam_data/Uc')
Uc_orig = Uc(:,1:100); 
% Bi
load('Beam_design/Bi_nom')
Ub_nom = Bi; 

load('Beam_design/Bi_opt')
Ub_opt = Bi; 

x_highfidelity = x_highfidelity';

% Plot Range to have indication: 
Uc_mean = mean(Uc_orig,2)';
Uc_max = max(Uc_orig,[],2)';
Uc_min = min(Uc_orig,[],2)';
Uc_yy = [Uc_min,fliplr(Uc_max)];

Uf_mean = mean(Uf,2)'; 
Uf_max = max(Uf,[],2)';
Uf_min = min(Uf,[],2)';
Uf_yy = [Uf_min,fliplr(Uf_max)];

Ucc_max = max(Uc_trial,[],2)';
Ucc_min = min(Uc_trial,[],2)';
Ucc_yy = [Ucc_min,fliplr(Ucc_max)];

Ub_nom_max = max(Ub_nom,[],2)';
Ub_nom_min = min(Ub_nom,[],2)';
Ub_nom_yy = [Ub_nom_min,fliplr(Ub_nom_max)];

Ub_opt_max = max(Ub_opt,[],2)';
Ub_opt_min = min(Ub_opt,[],2)';
Ub_opt_yy = [Ub_opt_min,fliplr(Ub_opt_max)];

xx = [x_highfidelity,fliplr(x_highfidelity)];


% Plot high, nominal low, and nominal bi ensembles
figure
hold on
h1 = fill(xx/50,Uf_yy,c1);
set(h1,'facealpha',.5)

h2 = fill(xx/50,Uc_yy,c2);
set(h2,'facealpha',.5)

h3 = fill(xx/50,Ub_nom_yy,c3);
set(h3,'facealpha',.5)

xlabel('x','interpreter','latex','Fontsize',FS)
ylabel('y','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([h1,h2,h3],{'H','Nominal L','Nominal B'},'interpreter', 'latex', 'fontsize', FS_leg)
ylim([-4, 0])
% % title('Gaussian')
grid on; 



% Plot high, nominal low, and nominal bi, optimal low and optimal bi ensembles
figure
hold on
h1 = fill(xx/50,Uf_yy,c1);
set(h1,'facealpha',.5)

h2 = fill(xx/50,Uc_yy,c2);
set(h2,'facealpha',.5)

h3 = fill(xx/50,Ub_nom_yy,c3);
set(h3,'facealpha',.5)

h4 = fill(xx/50,Ucc_yy,c4);
set(h4,'facealpha',.5)
p4 = plot(x_highfidelity/50,Ucc_min,'Color',c4,'LineWidth',2);

h5 = fill(xx/50,Ub_opt_yy,c5);
set(h5,'facealpha',.5)

xlabel('x','interpreter','latex','Fontsize',FS)
ylabel('y','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([h1,h2,h3,h4,h5],{'H','Nominal L','Nominal B','Optimal L','Optimal B'},'interpreter', 'latex', 'fontsize', FS_leg)
ylim([-4, 0])
% % title('Gaussian')
grid on; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot realization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify worst performing bi-fidelity estimate? 
% For now just the first sample

plot_index = 1; 

figure
hold on;
p1 = plot(x_highfidelity/50,Uf(:,plot_index),'Color',c1,'LineWidth', LW); 
p2 = plot(x_highfidelity/50,Uc_orig(:,plot_index),':','Color',c2,'LineWidth', LW); 
p3 = plot(x_highfidelity/50,Ub_nom(:,plot_index),'--','Color',c3,'LineWidth', LW); 
hold off
xlabel('x','interpreter','latex','Fontsize',FS)
ylabel('y','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1,p2,p3],{'H','Nominal L','Nominal B'},'interpreter', 'latex', 'fontsize', FS_leg)
ylim([-4, 0])
title('Displacement Realization' ,'interpreter', 'latex')
grid on; 

figure
hold on;
p1 = plot(x_highfidelity/50,Uf(:,plot_index),'Color',c1,'LineWidth', LW); 
p2 = plot(x_highfidelity/50,Uc_orig(:,plot_index),':','Color',c2,'LineWidth', LW); 
p3 = plot(x_highfidelity/50,Ub_nom(:,plot_index),'--','Color',c3,'LineWidth', LW); 
p4 = plot(x_highfidelity/50,Uc_trial(:,plot_index),':','Color',c4,'LineWidth', LW); 
p5 = plot(x_highfidelity/50,Ub_opt(:,plot_index),'-.','Color',c5,'LineWidth', LW); 
hold off
xlabel('x','interpreter','latex','Fontsize',FS)
ylabel('y','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1,p2,p3,p4,p5],{'H','Nominal L','Nominal B','Optimal L','Optimal B'},'interpreter', 'latex', 'fontsize', FS_leg)
ylim([-4, 0])
title('Displacement Realization' ,'interpreter', 'latex')
grid on; 

