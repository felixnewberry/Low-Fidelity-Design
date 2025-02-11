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

size_1 = [0,0,575,445]; 
size_2 = [0,0,1150,445]; 

% size_1 = [0,0,670,515]; 
% size_2 = [0,0,1340,515]; 

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
%%% Plot individual sensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% N = 2
load('Beam_design/individual_w_pm95')
error_w = error_bound_mat;
load('Beam_design/individual_h1_pm95')
error_h1 = error_bound_mat;
load('Beam_design/individual_h2_pm95')
error_h2 = error_bound_mat;
load('Beam_design/individual_h3_pm95')
error_h3 = error_bound_mat;

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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on

% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/beam_individual_all','epsc')
end

load('Beam_design/individual_h1h2_p2m1')
error_h1h2 = error_bound_mat;
delta_vec_h1h2 = delta_vec;

load('Beam_design/individual_h3_p35')
error_h3 = error_bound_mat;
delta_vec_h3 = delta_vec;

error_bound_nom = error_h3(1); 

figure
plot(100*delta_vec_h1h2,100*error_h1h2,'o-','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
xlabel('$ \Delta h_1 = \Delta h_2 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
ylim(100*[min([error_h3;error_h1h2]), max([error_h3;error_h1h2])])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_1)
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity $h_1 = h_2$','Interpreter','latex')


if save_on ==1
    saveas(gcf,'Plots/beam_individual_h1h2','epsc')
end

figure
plot(100*delta_vec_h3,100*error_h3,'o-','Color',c1,'LineWidth',LW,'MarkerSize',MS); 
xlabel('$ \Delta h_3 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
ylim(100*[min([error_h3;error_h1h2]), max([error_h3;error_h1h2])])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_1)
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity $h_3$','Interpreter','latex')

if save_on ==1
    saveas(gcf,'Plots/beam_individual_h3','epsc')
end


1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot grid search response surface 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Beam_design/grid_search')

% Find minimum bound
[min_bound,minIdx_b]=min(error_bound_mat(:));
[i_t1_B,i_t3_B]=ind2sub(size(error_bound_mat),minIdx_b);
t1_bound = delta_t1_vec(i_t1_B);
t3_bound = delta_t3_vec(i_t3_B);
min_bi_bound = error_Bi_mat(i_t1_B,i_t3_B);

% Find minimum bi
[min_bi,minIdx_bi]=min(error_Bi_mat(:));
[i_t1_bi,i_t3_bi]=ind2sub(size(error_Bi_mat),minIdx_bi);
t1_bi = delta_t1_vec(i_t1_bi);
t3_bi = delta_t3_vec(i_t3_bi);
min_bound_bi = error_bound_mat(i_t1_bi,i_t3_bi);

n_sim = 100; % loaded Uc and Ub data also needs to be adjusted if this is changed from 100


%%%%%%%%%%%%%%%%%%%
[xrow,ycol] = meshgrid(...
            linspace(100*delta_t3_vec(1),100*delta_t3_vec(end),21),...
            linspace(100*delta_t1_vec(1),100*delta_t1_vec(end),21)...
          );
[xq,yq] = meshgrid(...
            linspace(100*delta_t3_vec(1),100*delta_t3_vec(end),2100),...
            linspace(100*delta_t1_vec(1),100*delta_t1_vec(end),2100)...
          );

figure
hold on
error_bound_matq = interp2(xrow,ycol,100*error_bound_mat,xq,yq,'cubic');
surf(xq,yq,error_bound_matq)
max_z = 100*max(error_bound_matq(:)); 
% set(h, 'edgecolor','none');
shading interp
view(0,90)
colorbar
% 5, 0.2
p1 = plot3(0,0,max_z,'ro','MarkerSize',8,'linewidth',LW);
p2 = plot3(100*t3_bound,100*t1_bound,max_z,'rx','MarkerSize',8,'linewidth',LW);
p3 = plot3(100*t3_bi,100*t1_bi,max_z,'rs','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Nominal','Optimal', 'True Optimal'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthWest')
xlabel('$\Delta h_3$ [\%]','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1 = \Delta h_2$ [\%]','interpreter','latex','Fontsize',FS)
c =colorbar;
c.TickLabelInterpreter = 'latex'; 
axis tight
% caxis(100*[min([error_bound_mat(:); error_Bi_mat(:)]) max([error_bound_mat(:); error_Bi_mat(:)])])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)


%title('Error Bound','Interpreter','latex')

if save_on ==1
    saveas(gcf,'Plots/beam_grid_bound','epsc')
end


figure
hold on
error_Bi_matq = interp2(xrow,ycol,100*error_Bi_mat,xq,yq,'cubic');
surf(xq,yq,error_Bi_matq)
max_z = 100*max(error_Bi_matq(:)); 
shading interp
view(0,90)
colorbar
% 5, 0.2
p1 = plot3(0,0,max_z,'ro','MarkerSize',8,'linewidth',LW);
p2 = plot3(100*t3_bound,100*t1_bound,max_z,'rx','MarkerSize',8,'linewidth',LW);
p3 = plot3(100*t3_bi,100*t1_bi,max_z,'rs','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Nominal','Optimal', 'True Optimal'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthWest')
xlabel('$\Delta h_3$ [\%]','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1 = \Delta h_2$ [\%]','interpreter','latex','Fontsize',FS)
c =colorbar;
c.TickLabelInterpreter = 'latex'; 
axis tight
% caxis(100*[min([error_bound_mat(:); error_Bi_mat(:)]) max([error_bound_mat(:); error_Bi_mat(:)])])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/beam_grid_bi','epsc')
end

% Plot efficacy 
efficacy_mat = error_bound_mat./error_Bi_mat;

figure
hold on
efficacy_matq = interp2(xrow,ycol,100*efficacy_mat,xq,yq,'cubic');
surf(xq,yq,efficacy_matq)
max_z = 100*max(efficacy_matq(:));
shading interp
view(0,90)
colorbar
% 5, 0.2
p1 = plot3(0,0,max_z,'ro','MarkerSize',8,'linewidth',LW);
p2 = plot3(100*t3_bound,100*t1_bound,max_z,'rx','MarkerSize',8,'linewidth',LW);
p3 = plot3(100*t3_bi,100*t1_bi,max_z,'rs','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Nominal','Optimal', 'True Optimal'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthWest')
xlabel('$\Delta h_3$ [\%]','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1 = \Delta h_2$ [\%]','interpreter','latex','Fontsize',FS)
c =colorbar;
c.TickLabelInterpreter = 'latex'; 
axis tight
% caxis(100*[min([error_bound_mat(:); error_Bi_mat(:)]) max([error_bound_mat(:); error_Bi_mat(:)])])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)


if save_on ==1
    saveas(gcf,'Plots/beam_grid_efficacy','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot displacement ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Beam_data/x_highfidelity.txt')
load('Beam_data/xi')

load('Beam_data/Uf')
Uf = Uf(:,1:n_sim);

load('Beam_design/Nom')
Uc_nom = Uc;
Ub_nom = Ub; 

load('Beam_design/Opt')
Uc_opt = Uc;
Ub_opt = Ub; 

x_highfidelity = x_highfidelity';

% Plot Range to have indication: 
Uc_mean = mean(Uc_nom,2)';
Uc_max = max(Uc_nom,[],2)';
Uc_min = min(Uc_nom,[],2)';
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
h1 = fill(xx,Uf_yy,c1);
set(h1,'facealpha',.5)

h2 = fill(xx,Uc_yy,c2);
set(h2,'facealpha',.5)

h3 = fill(xx,Ub_nom_yy,c3);
set(h3,'facealpha',.5)

xlabel('$x$','interpreter','latex','Fontsize',FS)
ylabel('$y$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend([h1,h2,h3],{'HF','Nominal LF','Nominal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
ylim([-4, 0])
% % title('Gaussian')
grid on; 
set(gcf,'Position',size_1)
if save_on == 1
    saveas(gcf,'Plots/beam_ensemble_1','epsc')
end
% Plot high, nominal low, and nominal bi, optimal low and optimal bi ensembles
figure
hold on
h1 = fill(xx,Uf_yy,c1);
set(h1,'facealpha',.5)

h2 = fill(xx,Uc_yy,c2);
set(h2,'facealpha',.5)

h3 = fill(xx,Ub_nom_yy,c3);
set(h3,'facealpha',.5)

h4 = fill(xx,Ucc_yy,c4);
set(h4,'facealpha',.5)
p4 = plot(x_highfidelity/50,Ucc_min,'Color',c4,'LineWidth',2);

h5 = fill(xx,Ub_opt_yy,c5);
set(h5,'facealpha',.5)

xlabel('$x$','interpreter','latex','Fontsize',FS)
ylabel('$y$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend([h1,h2,h3,h4,h5],{'HF','Nominal LF','Nominal BF','Optimal LF','Optimal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
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

% Revert to standard 100
% Uf_all = Uf(:,1:3000);
Uf = Uf(:,1:n_sim);

%%% Standard 100 - as used for grid search and individual sensitivity
fprintf('Bound location and values: \n');
fprintf('t3: %d, t1=t2: %d \n',t3_bound, t1_bound);
fprintf('Bound: %d \n',min_bound);
fprintf('BF: %d \n',min_bi_bound);

fprintf('BF location and values \n');
fprintf('t3: %d, t1=t2: %d \n',t3_bi, t1_bi);
fprintf('Bound: %d \n',min_bound_bi);
fprintf('BF: %d \n',min_bi);

error_low_nom = norm(Uc_nom - Uf)/norm(Uf); 
error_low_opt = norm(Uc_opt - Uf)/norm(Uf); 
% error_bound_nom = norm(Ub_nom - Uf)/norm(Uf); % check was done with n = 2. 
error_bi_nom = norm(Ub_nom - Uf)/norm(Uf); % check was done with n = 2. 

% error_low_nom = 

fprintf('Results table nsim = 100: \n');
results = [error_low_nom, error_bound_nom, error_bi_nom; ...
    error_low_opt, min_bound, min_bi_bound]; 

results_tab = array2table(100*results, 'VariableNames',{'LF','Bound','BF'},...
    'RowNames',{'Nominal', 'Optimal'})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot displacement realization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify worst performing bi-fidelity estimate? 
% For now just the first sample

plot_index = 1; 

figure
hold on;
p1 = plot(x_highfidelity,Uf(:,plot_index),'Color',c1,'LineWidth', LW); 
p2 = plot(x_highfidelity,Uc_nom(:,plot_index),':','Color',c2,'LineWidth', LW); 
p3 = plot(x_highfidelity,Ub_nom(:,plot_index),'--','Color',c3,'LineWidth', LW); 
hold off
xlabel('$x$','interpreter','latex','Fontsize',FS)
ylabel('$y$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend([p1,p2,p3],{'HF','Nominal LF','Nominal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
ylim([-4, 0])
%title('Displacement Realization' ,'interpreter', 'latex')
% grid on; 
set(gcf,'Position',size_1)
if save_on ==1
    saveas(gcf,'Plots/beam_realization_1','epsc')
end

figure
hold on;
p1 = plot(x_highfidelity,Uf(:,plot_index),'Color',c1,'LineWidth', LW); 
p2 = plot(x_highfidelity,Uc_nom(:,plot_index),':','Color',c2,'LineWidth', LW); 
p3 = plot(x_highfidelity,Ub_nom(:,plot_index),'--','Color',c3,'LineWidth', LW); 
p4 = plot(x_highfidelity,Uc_opt(:,plot_index),':','Color',c4,'LineWidth', LW); 
p5 = plot(x_highfidelity,Ub_opt(:,plot_index),'-.','Color',c5,'LineWidth', LW); 
hold off
xlabel('$x$','interpreter','latex','Fontsize',FS)
ylabel('$y$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend([p1,p2,p3,p4,p5],{'HF','Nominal LF','Nominal BF','Optimal LF','Optimal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
ylim([-4, 0])
%title('Displacement Realization' ,'interpreter', 'latex')
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

load('Beam_design/tip_error_nom'); 
tip_error_nom_bi = tip_error_bi; 
tip_error_nom_L = tip_error_L; 

load('Beam_design/tip_error_opt'); 
tip_error_opt_bi = tip_error_bi; 
tip_error_opt_L = tip_error_L; 

figure
hold on
l1 = xline(error_bound_nom*100,'LineWidth', LW, 'Color', c1);
l2 = xline(min_bound*100,'--','LineWidth', LW, 'Color', c2);

h1 = histogram(abs(100*tip_error_nom_bi),n_hist,'FaceColor',c1);
h2 = histogram(abs(100*tip_error_opt_bi),n_hist,'FaceColor',c2);
hold off
legend([l1, l2, h1,h2],{'Nominal Bound','Optimal Bound','Nominal Ensemble','Optimal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
yl = ylim;
ylim([yl(1),yl(2)+1]);
xl = xlim;
xlim([xl(1),xl(2)+xl(2)*0.05]);

set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/beam_tip_bi_hist','epsc')
end

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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'Plots/beam_sv','epsc')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Sorted bound 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[bound_sort_s, bound_index_s] = sort(error_bound_mat(:));
bi_sort_s = error_Bi_mat(bound_index_s); 

delta_mat_s = delta_mat; 

figure
p1 = plot(100*bound_sort_s,'ob','color', c1, 'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(100*bi_sort_s,'xr', 'color', c2,'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot([1,200],[nom_bound_s_field*100, nom_bound_s_field*100],'b-','color', c1, 'LineWidth',LW); 
% p4 = plot([1,200],[nom_bi_s_field*100, nom_bi_s_field*100],'r--', 'color', c2,'LineWidth',LW); 
hold off
xlabel('Sample','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1,p2,p3,p4],{'Bound Samples','BF Samples','Nominal Bound','Nominal BF'},'interpreter', 'latex', 'fontsize', FS_leg, 'Location', 'NorthWest')
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
ylim([0,30])
% title('Stress Field','Interpreter','latex')
set(gcf,'Position',size_1)
if save_on ==1
    saveas(gcf,'plots/Beam_rand_samples_s','epsc')
end
