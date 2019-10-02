clear all
close all
clc

save_on = 1; 

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
%%% Plot line search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Beam_design/line_w_pm95')
error_w = error_bound_mat;
load('Beam_design/line_h1_pm95')
error_h1 = error_bound_mat;
load('Beam_design/line_h2_pm95')
error_h2 = error_bound_mat;
load('Beam_design/line_h3_pm95')
error_h3 = error_bound_mat;
load('Beam_design/delta_vec_pm95')

figure
hold on
p1 = plot(100*delta_vec,100*error_w,'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_vec,100*error_h1,'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_h2,'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_vec,100*error_h3,'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$\Delta [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1,p2,p3,p4],{'$w$','$h_1$','$h_2$','$h_3$'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','SouthEast')
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/line_all','epsc')
end

load('Beam_design/line_h1h2_p2m1')
error_h1h2 = error_bound_mat;
load('Beam_design/delta_vec_h1h2_p2m1')
delta_vec_h1h2 = delta_vec;

figure
plot(100*delta_vec_h1h2,100*error_h1h2,'o-','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
xlabel('$ \Delta h_1 = \Delta h_2 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/line_h1h2','epsc')
end

load('Beam_design/line_h3_p35')
error_h3 = error_bound_mat;
load('Beam_design/delta_vec_h3_p35')
delta_vec_h3 = delta_vec;

figure
plot(100*delta_vec_h3,100*error_h3,'o-','Color',c1,'LineWidth',LW,'MarkerSize',MS); 
xlabel('$ \Delta h_3 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
grid on
set(gcf,'Position',size_1)


if save_on ==1
    saveas(gcf,'Plots/line_h3','epsc')
end


figure 

subplot(1,2,1)
plot(100*delta_vec_h1h2,100*error_h1h2,'o-','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
xlabel('$ \Delta h_1 = \Delta h_2 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
ylim(100*[min(error_h3), max(error_h1h2)])
title('$h1$ and $h2$' ,'interpreter', 'latex')

subplot(1,2,2)
plot(100*delta_vec_h3,100*error_h3,'o-','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
xlabel('$ \Delta h_3 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
ylim(100*[min(error_h3), max(error_h1h2)])
title('$h3$' ,'interpreter', 'latex')

set(gcf,'Position',size_2)

if save_on ==1
    saveas(gcf,'Plots/line_h1h2h3','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot displacement ensemble
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
set(gcf,'Position',size_1)
saveas(gcf,'Plots/ensemble_1','epsc')

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
set(gcf,'Position',size_1)
if save_on ==1
    saveas(gcf,'Plots/ensemble_2','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot displacement realization
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
set(gcf,'Position',size_1)
if save_on ==1
    saveas(gcf,'Plots/realization_1','epsc')
end

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
set(gcf,'Position',size_1)
if save_on ==1
    saveas(gcf,'Plots/realization_2','epsc')
end



