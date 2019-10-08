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

size_1 = [0,0,445,345]; 
size_2 = [0,0,1340,515]; 

% size_1 = [0,0,670,515]; 
% size_2 = [0,0,1340,515]; 

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


%%% 1 rep 

%%% N = 2
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
% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/beam_line_all','epsc')
end

load('Beam_design/line_h1h2_p2m1')
error_h1h2 = error_bound_mat;
load('Beam_design/delta_vec_h1h2_p2m1')
delta_vec_h1h2 = delta_vec;

load('Beam_design/line_h3_p35')
error_h3 = error_bound_mat;
load('Beam_design/delta_vec_h3_p35')
delta_vec_h3 = delta_vec;

error_bound_nom = error_h3(1); 

figure
plot(100*delta_vec_h1h2,100*error_h1h2,'o-','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
xlabel('$ \Delta h_1 = \Delta h_2 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
ylim(100*[min([error_h3;error_h1h2]), max([error_h3;error_h1h2])])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
% grid on
set(gcf,'Position',size_1)
title('Line search $h_1 = h_2$','Interpreter','latex')

if save_on ==1
    saveas(gcf,'Plots/beam_line_h1h2','epsc')
end

figure
plot(100*delta_vec_h3,100*error_h3,'o-','Color',c1,'LineWidth',LW,'MarkerSize',MS); 
xlabel('$ \Delta h_3 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
ylim(100*[min([error_h3;error_h1h2]), max([error_h3;error_h1h2])])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
% grid on
set(gcf,'Position',size_1)
title('Line search $h_3$','Interpreter','latex')

if save_on ==1
    saveas(gcf,'Plots/beam_line_h3','epsc')
end


1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot grid search response surface 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('Beam_design/Grid_search_N10')
load('Beam_design/grid_search_N3')

figure
hold on
contourf(100*delta_t3_vec,100*delta_t1_vec,100*error_bound_mat)
colorbar
% 5, 0.2
p1 = plot(0,0,'ro','MarkerSize',8,'linewidth',LW);
p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',8,'linewidth',LW);
p3 = plot(100*t3_Bi,100*t1_Bi,'rs','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Nominal','Optimal', 'Bi'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthWest')
xlabel('$\Delta h_3$ [\%]','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1 = \Delta h_2$ [\%]','interpreter','latex','Fontsize',FS)
axis tight
caxis(100*[min([error_bound_mat(:); error_Bi_mat(:)]) max([error_bound_mat(:); error_Bi_mat(:)])])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
% colormap('cividis') % cividis, inferno, fake_parula, magma, plasma, viridis 
title('Error bound','Interpreter','latex')

if save_on ==1
    saveas(gcf,'Plots/beam_grid_bound','epsc')
end

figure
hold on
contourf(100*delta_t3_vec,100*delta_t1_vec,100*error_Bi_mat)
colorbar
p1 = plot(0,0,'ro','MarkerSize',8,'linewidth',LW);
p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',8,'linewidth',LW);
p3 = plot(100*t3_Bi,100*t1_Bi,'rs','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Nominal','Optimal', 'Bi'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthWest')
xlabel('$\Delta h_3 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1= \Delta h_2$ [\%]','interpreter','latex','Fontsize',FS)
axis tight
caxis(100*[min([error_bound_mat(:); error_Bi_mat(:)]) max([error_bound_mat(:); error_Bi_mat(:)])])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
% colormap('cividis')
title('Bi-fidelity error','Interpreter','latex')

if save_on ==1
    saveas(gcf,'Plots/beam_grid_bi','epsc')
end


% data

% A_hat
% Minimum of each column
[min_values, index_t1_A] = min(error_Bi_mat); 
% Global minimum, column index is t3
[min_value, index_t3_A] = min(min_values); 
index_t1_A = index_t1_A(index_t3_A);
t1_Bi = delta_t1_vec(index_t1_A);
t3_Bi = delta_t3_vec(index_t3_A);
min_Bi = min_value; 

min_bound_Bi = error_bound_mat(index_t1_A,index_t3_A);

% Error bound
% Minimum of each column
[min_values, index_t1_B] = min(error_bound_mat); 
% Global minimum, column index is t3
[min_value, index_t3_B] = min(min_values); 
index_t1_B = index_t1_B(index_t3_B);
t1_bound = delta_t1_vec(index_t1_B);
t3_bound = delta_t3_vec(index_t3_B);
min_bound = min_value; 

min_Bi_bound = error_Bi_mat(index_t1_B,index_t3_B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot displacement ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Beam_data/x_highfidelity.txt')
load('Beam_data/xi')

load('Beam_design/Uc_opt')
Uc_opt = Uc;

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

Ucc_max = max(Uc_opt,[],2)';
Ucc_min = min(Uc_opt,[],2)';
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
saveas(gcf,'Plots/beam_ensemble_1','epsc')

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
    saveas(gcf,'Plots/beam_ensemble_2','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Results table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Bound location and values: \n');
fprintf('t3: %d, t1=t2: %d \n',t3_bound, t1_bound);
fprintf('Bound: %d \n',min_bound);
fprintf('Bi: %d \n',min_Bi_bound);

fprintf('Bi location and values \n');
fprintf('t3: %d, t1=t2: %d \n',t3_Bi, t1_Bi);
fprintf('Bound: %d \n',min_bound_Bi);
fprintf('Bi: %d \n',min_Bi);

1; 

error_low_nom = norm(Uc_orig - Uf)/norm(Uf); 
error_low_opt = norm(Uc_opt - Uf)/norm(Uf); 
% error_bound_nom = norm(Ub_nom - Uf)/norm(Uf); % check was done with n = 2. 
error_bi_nom = norm(Ub_nom - Uf)/norm(Uf); % check was done with n = 2. 

% error_low_nom = 

fprintf('Results table: \n');
results = [error_low_nom, error_bound_nom, error_bi_nom; ...
    error_low_opt, min_bound, min_Bi_bound]; 

results_tab = array2table(100*results, 'VariableNames',{'Low','Bound','Bi'},...
    'RowNames',{'Nominal', 'Optimal'})

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
% grid on; 
set(gcf,'Position',size_1)
if save_on ==1
    saveas(gcf,'Plots/beam_realization_1','epsc')
end

figure
hold on;
p1 = plot(x_highfidelity/50,Uf(:,plot_index),'Color',c1,'LineWidth', LW); 
p2 = plot(x_highfidelity/50,Uc_orig(:,plot_index),':','Color',c2,'LineWidth', LW); 
p3 = plot(x_highfidelity/50,Ub_nom(:,plot_index),'--','Color',c3,'LineWidth', LW); 
p4 = plot(x_highfidelity/50,Uc_opt(:,plot_index),':','Color',c4,'LineWidth', LW); 
p5 = plot(x_highfidelity/50,Ub_opt(:,plot_index),'-.','Color',c5,'LineWidth', LW); 
hold off
xlabel('x','interpreter','latex','Fontsize',FS)
ylabel('y','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1,p2,p3,p4,p5],{'H','Nominal L','Nominal B','Optimal L','Optimal B'},'interpreter', 'latex', 'fontsize', FS_leg)
ylim([-4, 0])
title('Displacement Realization' ,'interpreter', 'latex')
% grid on; 
set(gcf,'Position',size_1)
if save_on ==1
    saveas(gcf,'Plots/beam_realization_2','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot tip displacement histograms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Vizualize  data
n_hist = 20; 

% I think some of the samples suck. look at 647 on... 

load('Beam_design/tip_error_nom_bi'); 
tip_error_nom_bi = tip_error; 
load('Beam_design/tip_error_nom_L')
tip_error_nom_L = tip_error_L; 

load('Beam_design/tip_error_opt_bi'); 
tip_error_opt_bi = tip_error; 
load('Beam_design/tip_error_opt_L')
tip_error_opt_L = tip_error_L; 

figure
hold on
h1 = histogram(abs(tip_error_nom_bi),n_hist,'FaceColor',c1);
h2 = histogram(abs(tip_error_opt_bi),n_hist,'FaceColor',c2);
hold off
legend([h1,h2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('bi-fidelity error','interpreter','latex','Fontsize',FS)
ylabel('frequency','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/beam_tip_bi_hist','epsc')
end
%%% Optimal L is a very narrow distribution - looks silly. 
% figure
% hold on
% h1 = histogram(abs(tip_error_nom_L),n_hist);
% h2 = histogram(abs(tip_error_opt_L),n_hist);
% hold off
% legend([h1,h2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
% xlabel('bi-fidelity error','interpreter','latex','Fontsize',FS)
% ylabel('frequency','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot eigenvalue decay 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1; 

sv_L_nom = svd(Uc); 
sv_L_Opt = svd(Uc_opt); 

i_lim = 5; 

% normalized
figure
p1 = semilogy(sv_L_nom(1:i_lim)/sv_L_nom(1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = semilogy(sv_L_Opt(1:i_lim)/sv_L_Opt(1),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('Index $i$','interpreter','latex','Fontsize',FS)
ylabel('Normalized Singular Values','interpreter','latex','Fontsize',FS)
legend([p1,p2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/beam_sv_normalized','epsc')
end

% normalized by nominal 
figure
p1 = semilogy(sv_L_nom(1:i_lim)/sv_L_nom(1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = semilogy(sv_L_Opt(1:i_lim)/sv_L_nom(1),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('Index $i$','interpreter','latex','Fontsize',FS)
ylabel('Singular Values','interpreter','latex','Fontsize',FS)
legend([p1,p2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/beam_sv_nom','epsc')
end

% not normalized 
figure
p1 = semilogy(sv_L_nom(1:i_lim),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = semilogy(sv_L_Opt(1:i_lim),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('Index $i$','interpreter','latex','Fontsize',FS)
ylabel('Singular Values','interpreter','latex','Fontsize',FS)
legend([p1,p2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/beam_sv','epsc')
end
