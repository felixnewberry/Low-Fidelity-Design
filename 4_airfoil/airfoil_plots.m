% Airfoil Plots

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

% Paper size
size_1 = [0,0,575,445]; 
size_2 = [0,0,1150,445]; 

size_square = [0,0,445,445]; 
size_large = [0,0,668,518]; 


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
%%% r calibration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('airfoil_data/Uf.mat', 'Uf')
load('airfoil_data/Uc_nom.mat', 'Uc_nom')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Individual Sensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load('results_plot/plot_data_pm50_r10')
% load('results_plot/plot_data_pm10_r10')

load('results_plot/plot_data_pm50_r6')
load('results_plot/plot_data_pm10_r6')

bound_nom = bound_plot_pm50(1,11); bi_nom = bi_plot_pm50(1,11); low_nom = err_low_plot_pm50(1,11); 


figure
p1 = plot(delta_vec_pm10, 100*bound_plot_pm10(1,:), 'o-','LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(delta_vec_pm10, 100*bound_plot_pm10(2,:), 'x:','LineWidth',LW,'MarkerSize',MS); 
p3 = plot(delta_vec_pm10, 100*bound_plot_pm10(3,:), '+-.','LineWidth',LW,'MarkerSize',MS); 
p4 = plot(delta_vec_pm10, 100*bound_plot_pm10(4,:), '*--','LineWidth',LW,'MarkerSize',MS); 
p5 = plot(delta_vec_pm10, 100*bound_plot_pm10(5,:), 's-','LineWidth',LW,'MarkerSize',MS); 
p6 = plot(delta_vec_pm10, 100*bound_plot_pm10(6,:), 'd:','LineWidth',LW,'MarkerSize',MS); 
p7 = plot(delta_vec_pm10, 100*bound_plot_pm10(7,:), 'v-.','LineWidth',LW,'MarkerSize',MS); 
p8 = plot(delta_vec_pm10, 100*bound_plot_pm10(8,:), '^--','LineWidth',LW,'MarkerSize',MS); 
p9 = plot(delta_vec_pm10, 100*bound_plot_pm10(9,:), '<-','LineWidth',LW,'MarkerSize',MS); 
p10 = plot(delta_vec_pm10, 100*bound_plot_pm10(10,:), '>:','LineWidth',LW,'MarkerSize',MS); 
p11 = plot(delta_vec_pm10, 100*bound_plot_pm10(11,:), 'p-.','LineWidth',LW,'MarkerSize',MS); 
p12 = plot(delta_vec_pm10, 100*bound_plot_pm10(12,:), 'h--','LineWidth',LW,'MarkerSize',MS); 
legend([p1 ,p2, p3, p4 ,p5, p6, p7 ,p8, p9, p10 ,p11, p12],{'$\rho$','$m$', '$p$', '$t$', '$\alpha$','$c_{b1}$','$c_{b2}$', '$c_{w2}$', '$c_{w3}$', '$\kappa$', '$\sigma$', '$c_{v1}$'},'interpreter', 'latex', 'fontsize', FS_leg)

xlabel('$\Delta[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'plots/airfoil_line_pm10_r6','epsc')
end

figure
subplot(1,2,1)
p1 = plot(delta_vec_pm10, 100*bound_plot_pm10(1,:), 'o-','LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(delta_vec_pm10, 100*bound_plot_pm10(2,:), 'x:','LineWidth',LW,'MarkerSize',MS); 
p3 = plot(delta_vec_pm10, 100*bound_plot_pm10(3,:), '+-.','LineWidth',LW,'MarkerSize',MS); 
p4 = plot(delta_vec_pm10, 100*bound_plot_pm10(4,:), '*--','LineWidth',LW,'MarkerSize',MS); 
p5 = plot(delta_vec_pm10, 100*bound_plot_pm10(5,:), 's-','LineWidth',LW,'MarkerSize',MS); 
p6 = plot(delta_vec_pm10, 100*bound_plot_pm10(6,:), 'd:','LineWidth',LW,'MarkerSize',MS); 
p7 = plot(delta_vec_pm10, 100*bound_plot_pm10(7,:), 'v-.','LineWidth',LW,'MarkerSize',MS); 
p8 = plot(delta_vec_pm10, 100*bound_plot_pm10(8,:), '^--','LineWidth',LW,'MarkerSize',MS); 
p9 = plot(delta_vec_pm10, 100*bound_plot_pm10(9,:), '<-','LineWidth',LW,'MarkerSize',MS); 
p10 = plot(delta_vec_pm10, 100*bound_plot_pm10(10,:), '>:','LineWidth',LW,'MarkerSize',MS); 
p11 = plot(delta_vec_pm10, 100*bound_plot_pm10(11,:), 'p-.','LineWidth',LW,'MarkerSize',MS); 
p12 = plot(delta_vec_pm10, 100*bound_plot_pm10(12,:), 'h--','LineWidth',LW,'MarkerSize',MS); 
legend([p1 ,p2, p3, p4 ,p5, p6, p7 ,p8, p9, p10 ,p11, p12],{'$\rho$','$m$', '$p$', '$t$', '$\alpha$','$c_{b1}$','$c_{b2}$', '$c_{w2}$', '$c_{w3}$', '$\kappa$', '$\sigma$', '$c_{v1}$'},'interpreter', 'latex', 'fontsize', FS_leg)

xlabel('$\Delta[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on


subplot(1,2,2)
p1 = plot(delta_vec_pm10, 100*bi_plot_pm10(1,:), 'o-','LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(delta_vec_pm10, 100*bi_plot_pm10(2,:), 'x:','LineWidth',LW,'MarkerSize',MS); 
p3 = plot(delta_vec_pm10, 100*bi_plot_pm10(3,:), '+-.','LineWidth',LW,'MarkerSize',MS); 
p4 = plot(delta_vec_pm10, 100*bi_plot_pm10(4,:), '*--','LineWidth',LW,'MarkerSize',MS); 
p5 = plot(delta_vec_pm10, 100*bi_plot_pm10(5,:), 's-','LineWidth',LW,'MarkerSize',MS); 
p6 = plot(delta_vec_pm10, 100*bi_plot_pm10(6,:), 'd:','LineWidth',LW,'MarkerSize',MS); 
p7 = plot(delta_vec_pm10, 100*bi_plot_pm10(7,:), 'v-.','LineWidth',LW,'MarkerSize',MS); 
p8 = plot(delta_vec_pm10, 100*bi_plot_pm10(8,:), '^--','LineWidth',LW,'MarkerSize',MS); 
p9 = plot(delta_vec_pm10, 100*bi_plot_pm10(9,:), '<-','LineWidth',LW,'MarkerSize',MS); 
p10 = plot(delta_vec_pm10, 100*bi_plot_pm10(10,:), '>:','LineWidth',LW,'MarkerSize',MS); 
p11 = plot(delta_vec_pm10, 100*bi_plot_pm10(11,:), 'p-.','LineWidth',LW,'MarkerSize',MS); 
p12 = plot(delta_vec_pm10, 100*bi_plot_pm10(12,:), 'h--','LineWidth',LW,'MarkerSize',MS); 

xlabel('$\Delta[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Bi Error $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1 ,p2, p3, p4 ,p5, p6, p7 ,p8, p9, p10 ,p11, p12],{'$\rho$','$m$', '$p$', '$t$', '$\alpha$','$c_{b1}$','$c_{b2}$', '$c_{w2}$', '$c_{w3}$', '$\kappa$', '$\sigma$', '$c_{v1}$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_2)


figure
p1 = plot(delta_vec_pm50, 100*bound_plot_pm50(1,:), 'o-','LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(delta_vec_pm50, 100*bound_plot_pm50(2,:), 'x:','LineWidth',LW,'MarkerSize',MS); 
p3 = plot(delta_vec_pm50, 100*bound_plot_pm50(3,:), '+-.','LineWidth',LW,'MarkerSize',MS); 
p4 = plot(delta_vec_pm50, 100*bound_plot_pm50(4,:), '*--','LineWidth',LW,'MarkerSize',MS); 
p5 = plot(delta_vec_pm50, 100*bound_plot_pm50(5,:), 's-','LineWidth',LW,'MarkerSize',MS); 
p6 = plot(delta_vec_pm50, 100*bound_plot_pm50(6,:), 'd:','LineWidth',LW,'MarkerSize',MS); 
legend([p1 ,p2, p3, p4 ,p5, p6],{'$\rho$','$m$', '$p$', '$t$', '$\alpha$','$c_{b1}$'},'interpreter', 'latex', 'fontsize', FS_leg)

xlabel('$\Delta[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'plots/airfoil_line_pm50_r6','epsc')
end
figure
subplot(1,2,1)
p1 = plot(delta_vec_pm50, 100*bound_plot_pm50(1,:), 'o-','LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(delta_vec_pm50, 100*bound_plot_pm50(2,:), 'x:','LineWidth',LW,'MarkerSize',MS); 
p3 = plot(delta_vec_pm50, 100*bound_plot_pm50(3,:), '+-.','LineWidth',LW,'MarkerSize',MS); 
p4 = plot(delta_vec_pm50, 100*bound_plot_pm50(4,:), '*--','LineWidth',LW,'MarkerSize',MS); 
p5 = plot(delta_vec_pm50, 100*bound_plot_pm50(5,:), 's-','LineWidth',LW,'MarkerSize',MS); 
p6 = plot(delta_vec_pm50, 100*bound_plot_pm50(6,:), 'd:','LineWidth',LW,'MarkerSize',MS); 
legend([p1 ,p2, p3, p4 ,p5, p6],{'$\rho$','$m$', '$p$', '$t$', '$\alpha$','$c_{b1}$'},'interpreter', 'latex', 'fontsize', FS_leg)

xlabel('$\Delta[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on


subplot(1,2,2)
p1 = plot(delta_vec_pm50, 100*bi_plot_pm50(1,:), 'o-','LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(delta_vec_pm50, 100*bi_plot_pm50(2,:), 'x:','LineWidth',LW,'MarkerSize',MS); 
p3 = plot(delta_vec_pm50, 100*bi_plot_pm50(3,:), '+-.','LineWidth',LW,'MarkerSize',MS); 
p4 = plot(delta_vec_pm50, 100*bi_plot_pm50(4,:), '*--','LineWidth',LW,'MarkerSize',MS); 
p5 = plot(delta_vec_pm50, 100*bi_plot_pm50(5,:), 's-','LineWidth',LW,'MarkerSize',MS); 
p6 = plot(delta_vec_pm50, 100*bi_plot_pm50(6,:), 'd:','LineWidth',LW,'MarkerSize',MS); 

xlabel('$\Delta[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Bi Error $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1 ,p2, p3, p4 ,p5, p6],{'$\rho$','$m$', '$p$', '$t$', '$\alpha$','$c_{b1}$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 500 samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


design_params_500 = importdata(strcat('results_data', '/design_params_500.dat'), ',' , 1); 
design_params_500 = design_params_500.data; 

% load('batch_500_res_r10', 'cp_results_500', 'err_bound_vec_500', 'err_bi_vec_500','err_low_vec_500', 'N', 'n', 'r'); 
load('cp_results_500')
cp_results_500 = cp_results; 
load('batch_500_res', 'err_bound_vec_500', 'err_bi_vec_500','err_low_vec_500', 'N', 'n', 'r'); 

[bound_order, bound_i] = sort(err_bound_vec_500); 
[bi_order, bi_i] = sort(err_bi_vec_500); 

figure
p3 = plot(100*err_bound_vec_500(bound_i), 'o','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
hold on
p4 = plot(100*err_bi_vec_500(bound_i), 'x','Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p1 = plot(ones(500,1)*100*bound_nom, '-', 'Color',c3,'LineWidth',LW,'MarkerSize',MS); 

p2 = plot(ones(500,1)*100*bi_nom, '--','Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
legend([p3,p4, p1, p2],{'Bound Samples','BF Samples', 'Nominal Bound', 'Nominal BF'},'interpreter', 'latex', 'fontsize', FS_leg, 'Location', 'NorthWest')
hold off
xlabel('Sample','interpreter','latex','Fontsize',FS)
ylabel('Error $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on

if save_on ==1
    saveas(gcf,'plots/airfoil_rand_samples_500_r6','epsc')
end

%%% Closer look - best 100 samples

%%% Closer look - best 100 samples
figure
p3 = plot(100*err_bound_vec_500(bound_i(1:100)), 'o','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
hold on
p4 = plot(100*err_bi_vec_500(bound_i(1:100)), 'x','Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p1 = plot(100*ones(100,1)*bound_nom, '-', 'Color',c3,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*ones(100,1)*bi_nom, '--','Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
legend([p3,p4, p1, p2],{'Bound Samples','BF Samples', 'Nominal Bound', 'Nominal BF'},'interpreter', 'latex', 'fontsize', FS_leg, 'Location', 'NorthWest')
hold off
xlabel('Sample','interpreter','latex','Fontsize',FS)
ylabel('Error $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
yl = ylim;
ylim([yl(1),yl(2)*(1+0.05)]);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on

res_mat_500 = 100*[low_nom, bound_nom, bi_nom; ...
    err_low_vec_500(bound_i(1)), err_bound_vec_500(bound_i(1)), err_bi_vec_500(bound_i(1)); ...
    err_low_vec_500(bi_i(1)), err_bound_vec_500(bi_i(1)), err_bi_vec_500(bi_i(1))];

res_tab_500 = array2table(res_mat_500, 'VariableNames', {'Low', 'Bound', 'Bi'}, 'RowNames', {'Nominal', 'Opt Bound', 'Opt Bi'}); 
res_tab_500

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Histogram 500
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cp_opt_bound_l = cp_results_500(:,:,bound_i(1)); 
cp_opt_bi_l = cp_results_500(:,:,bi_i(1)); 
cp_bound_nom = Uc_nom;
% I could check this data with where bound_i(1) came from...

%%%
% BF estimate 
rand_sample = 1:n; 

B_bound = cp_opt_bound_l/norm(cp_opt_bound_l,'fro');
A = Uf/norm(Uf,'fro');
[P_s,ix] = matrixIDvR(B_bound,r);
cp_bi_bound = A(:,ix)*P_s*norm(Uf,'fro'); 

B_bi = cp_opt_bi_l/norm(cp_opt_bi_l,'fro');
[P_s,ix] = matrixIDvR(B_bi,r);
cp_bi_bi = A(:,ix)*P_s*norm(Uf,'fro'); 

B_bi = cp_bound_nom/norm(cp_bound_nom,'fro');
[P_s,ix] = matrixIDvR(B_bi,r);
cp_bi_nom = A(:,ix)*P_s*norm(Uf,'fro');  
err_Bhat = norm(A-A(:,ix)*P_s) 

%%%%%%%%%%%% 

%%% Vizualize  data
n_hist = 20; 

error_b_nom = vecnorm(cp_bi_nom - Uf)./vecnorm(Uf);
error_b_opt = vecnorm(cp_bi_bound - Uf)./vecnorm(Uf);
error_b_bi = vecnorm(cp_bi_bi - Uf)./vecnorm(Uf);


% error_b_bi % looks very good... 
% error_b_opt %does not. 
% Check B_bi and B_bound

% check data on this one... 
figure
hold on
l1 = xline(res_mat_500(1,2),'LineWidth', LW, 'Color', c1);
l2 = xline(res_mat_500(2,2),'--','LineWidth', LW, 'Color', c2);

h1 = histogram(abs(100*error_b_nom),n_hist,'FaceColor',c1);
h2 = histogram(abs(100*error_b_opt),n_hist,'FaceColor',c2);
hold off
legend([l1, l2, h1,h2],{'Nominal Bound','Optimal Bound','Nominal Ensemble','Optimal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
yl = ylim;
ylim([yl(1),yl(2)*(1+0.05)]);
xl = xlim;
xlim([xl(1),xl(2)*(1+0.05)]);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'plots/airfoil_cp_hist_opt_bound_500_r6','epsc')
end

figure
hold on
l1 = xline(res_mat_500(1,2),'LineWidth', LW, 'Color', c1);
l2 = xline(res_mat_500(2,2),'--','LineWidth', LW, 'Color', c2);

h1 = histogram(abs(100*error_b_nom),n_hist,'FaceColor',c1);
h2 = histogram(abs(100*error_b_bi),n_hist,'FaceColor',c2);
hold off
legend([l1, l2, h1,h2],{'Nominal Bound','Optimal Bound','Nominal Ensemble','Optimal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
yl = ylim;
ylim([yl(1),yl(2)*(1+0.05)]);
xl = xlim;
xlim([xl(1),xl(2)*(1+0.05)]);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'plots/airfoil_cp_hist_opt_bi_500_r6','epsc')
end