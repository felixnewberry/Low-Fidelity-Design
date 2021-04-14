% Airfoil Plots

clear all
close all
clc
 
save_on = 0; 

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

load('results_data/cp_results_lims_pm5.mat', 'cp_results')
cp_results_pm5 = cp_results; 
design_params_pm5 = importdata(strcat('results_data', '/design_params_pm5.dat'), ',' , 1); 
design_params_pm5 = design_params_pm5.data; 

% pm10
load('results_data/cp_results_lims_pm10.mat', 'cp_results')
cp_results_pm10 = cp_results; 
design_params_pm10 = importdata(strcat('results_data', '/design_params_pm10.dat'), ',' , 1); 
design_params_pm10 = design_params_pm10.data; 

load('results_data/cp_results_lims_pm15_20.mat', 'cp_results')
cp_results_pm15_20 = cp_results; 
design_params_pm15_20 = importdata(strcat('results_data', '/design_params_pm15_20.dat'), ',' , 1); 
design_params_pm15_20 = design_params_pm15_20.data; 


load('results_data/cp_results_lims_pm25_30.mat', 'cp_results')
cp_results_pm25_30 = cp_results; 
design_params_pm25_30 = importdata(strcat('results_data', '/design_params_pm25_30.dat'), ',' , 1); 
design_params_pm25_30 = design_params_pm25_30.data; 


load('results_data/cp_results_lims_pm50.mat', 'cp_results')
cp_results_pm50 = cp_results; 
design_params_pm50 = importdata(strcat('results_data', '/design_params_pm50.dat'), ',' , 1); 
design_params_pm50 = design_params_pm50.data; 
% order is pm 0.35, 40, 45, 50 then through each param for first 6. 

% high-fidelity
load('airfoil_data/Uf.mat', 'Uf')

Uc_nom = cp_results_pm5(:,:,1); 

N = 500; 

%%% Calibrate

% Calibrate nominal bound 
r_vec = 3:30;
err_bound_vec_pm5 = zeros(1,length(r_vec)); 
err_bi_vec_pm5 = zeros(1,length(r_vec)); 

for i_r = 1:length(r_vec)
    [err_bound_vec_pm5(i_r),err_bi_vec_pm5(i_r), err_low] = my_airfoil_bound(Uc_nom,Uf, N, r_vec(i_r)+10, r_vec(i_r));
end

figure
p1 = plot(r_vec,100*err_bi_vec_pm5,'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(r_vec,100*err_bound_vec_pm5,'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$r$','interpreter','latex','Fontsize',FS)
ylabel('Error $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1,p2],{'BF','Bound'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on ==1
%     saveas(gcf,'plots_airfoil/airfoil_r_calibration','epsc')
    saveas(gcf,'plots_airfoil/airfoil_r_calibration','png')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Individual Sensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Line
% 
r = 6; 
n = r+10; 

% r = 10; 
% n = r+10; 

n_ens_pm5 = 25; 

err_bound_vec_pm5 = zeros(1,n_ens_pm5); 
err_bi_vec_pm5 = zeros(1,n_ens_pm5); 
err_low_vec_pm5 = zeros(1,n_ens_pm5); 
for i_ens = 1:n_ens_pm5
    [err_bound_vec_pm5(i_ens),err_bi_vec_pm5(i_ens), err_low_vec_pm5(i_ens)] = my_airfoil_bound(cp_results_pm5(:,:,i_ens),Uf, N, n, r);
end

n_ens_pm10 = 24; 

err_bound_vec_pm10 = zeros(1,n_ens_pm10); 
err_bi_vec_pm10 = zeros(1,n_ens_pm10); 
err_low_vec_pm10 = zeros(1,n_ens_pm10); 
for i_ens = 1:n_ens_pm10
    [err_bound_vec_pm10(i_ens),err_bi_vec_pm10(i_ens), err_low_vec_pm10(i_ens)] = my_airfoil_bound(cp_results_pm10(:,:,i_ens),Uf, N, n, r);
end

n_ens_pm15_20 = 24; 

err_bound_vec_pm15_20 = zeros(1,n_ens_pm15_20); 
err_bi_vec_pm15_20 = zeros(1,n_ens_pm15_20); 
err_low_vec_pm15_20 = zeros(1,n_ens_pm15_20); 

err_bound_vec_pm30 = zeros(1,n_ens_pm15_20); 
err_bi_vec_pm30 = zeros(1,n_ens_pm15_20); 
err_low_vec_pm30 = zeros(1,n_ens_pm15_20); 

for i_ens = 1:n_ens_pm15_20
    [err_bound_vec_pm15_20(i_ens),err_bi_vec_pm15_20(i_ens), err_low_vec_pm15_20(i_ens)] = my_airfoil_bound(cp_results_pm15_20(:,:,i_ens),Uf, N, n, r);
    [err_bound_vec_pm30(i_ens),err_bi_vec_pm30(i_ens), err_low_vec_pm30(i_ens)] = my_airfoil_bound(cp_results_pm25_30(:,:,i_ens),Uf, N, n, r);

end

n_ens_pm50 = 48; 

err_bound_vec_pm50 = zeros(1,n_ens_pm50); 
err_bi_vec_pm50 = zeros(1,n_ens_pm50); 
err_low_vec_pm50 = zeros(1,n_ens_pm50); 

for i_ens = 1:n_ens_pm50
    [err_bound_vec_pm50(i_ens),err_bi_vec_pm50(i_ens), err_low_vec_pm50(i_ens)] = my_airfoil_bound(cp_results_pm50(:,:,i_ens),Uf, N, n, r);
end


% 1 is nominal 
bound_nom = err_bound_vec_pm5(1); bi_nom = err_bi_vec_pm5(1); low_nom = err_low_vec_pm5(1); 

% density, m, p, t, a, sa1 - 7
i_4 = [2:2:25]; 
i_2 = [3:2:25]; 
i_1 = [2:2:24];
i_5 = [1:2:24];



bound_plot_pm10 = [err_bound_vec_pm10(i_1)', err_bound_vec_pm5(i_2)', ones(length(i_2),1)*bound_nom, err_bound_vec_pm5(i_4)', err_bound_vec_pm10(i_5)'];
bi_plot_pm10 = [err_bi_vec_pm10(i_1)', err_bi_vec_pm5(i_2)', ones(length(i_2),1)*bi_nom, err_bi_vec_pm5(i_4)', err_bi_vec_pm10(i_5)'];

% 
delta_vec_pm10 = [-10, -5, 0, 5, 10];

% save('results_plot/plot_data_pm10_r6', 'bound_plot_pm10', 'bi_plot_pm10', 'delta_vec_pm10');


%%% +- 10 
figure
subplot(1,2,1)
p2 = plot([-10, -5, 0, 5, 10], 100*bound_plot_pm10(1:5,:)', 'o-','LineWidth',LW,'MarkerSize',MS); 
hold on
p1 = plot([-10, -5, 0, 5, 10], 100*bound_plot_pm10(6:end,:)', 'x--','LineWidth',LW,'MarkerSize',MS); 
xlabel('Delta $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on



subplot(1,2,2)
% p1 = plot([-5, 0, 5], bound_plot(1:5,:), 'o--','LineWidth',LW,'MarkerSize',MS); 
% hold on
p1 = plot([-10, -5, 0, 5, 10], 100*bi_plot_pm10(1:5,:)', 'o-','LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot([-10, -5, 0, 5, 10], 100*bi_plot_pm10(6:end,:), 'x--','LineWidth',LW,'MarkerSize',MS); 

xlabel('Delta $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('BF Error $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1(1),p1(2), p1(3), p1(4), p1(5), p2(1), p2(2), p2(3), p2(4), p2(5), p2(6), p2(7)],{'$\rho$','m', 'p', 't', 'a','SA1','SA2', 'SA3', 'SA4', 'SA5', 'SA6', 'SA7'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_2)

if save_on ==1
    saveas(gcf,'plots_airfoil/airfoil_line_pm10_r10','epsc')
    saveas(gcf,'plots_airfoil/airfoil_line_pm10_r10','png')
end

%%% +- 50 

% Drop the 6 other SA parameters
i_6 = [2:2:13]; 
i_4 = [3:2:13]; 
i_3 = [2:2:12];
i_7 = [1:2:12]; 

i_1 = [14:2:24]; 
i_2 = [2:2:12]; 
i_8 = [1:2:12]; 
i_9 = [13:2:24]; 


i_35n = [2:8:48];
i_35p = [1:8:48];
i_40n = [4:8:48];
i_40p = [3:8:48];
i_45n = [6:8:48];
i_45p = [5:8:48];
i_50n = [8:8:48];
i_50p = [7:8:48];

bound_plot_pm50 = [err_bound_vec_pm50(i_50n)', err_bound_vec_pm50(i_45n)', ...
    err_bound_vec_pm50(i_40n)', err_bound_vec_pm50(i_35n)', err_bound_vec_pm30(i_1)', err_bound_vec_pm30(i_2)', err_bound_vec_pm15_20(i_1)', err_bound_vec_pm15_20(i_2)',...
    err_bound_vec_pm10(i_3)', err_bound_vec_pm5(i_4)', ones(length(i_4),1)*bound_nom, ...
    err_bound_vec_pm5(i_6)', err_bound_vec_pm10(i_7)', err_bound_vec_pm15_20(i_8)',...
    err_bound_vec_pm15_20(i_9)', err_bound_vec_pm30(i_8)', err_bound_vec_pm30(i_9)', ...
    err_bound_vec_pm50(i_35p)', err_bound_vec_pm50(i_40p)', err_bound_vec_pm50(i_45p)', err_bound_vec_pm50(i_50p)'];
bi_plot_pm50 = [err_bi_vec_pm50(i_50n)', err_bi_vec_pm50(i_45n)', err_bi_vec_pm50(i_40n)', ...
    err_bi_vec_pm50(i_35n)', err_bi_vec_pm30(i_1)', err_bi_vec_pm30(i_2)', err_bi_vec_pm15_20(i_1)', err_bi_vec_pm15_20(i_2)', err_bi_vec_pm10(i_3)',...
    err_bi_vec_pm5(i_4)', ones(length(i_4),1)*bi_nom, err_bi_vec_pm5(i_6)', ...
    err_bi_vec_pm10(i_7)', err_bi_vec_pm15_20(i_8)', err_bi_vec_pm15_20(i_9)',...
    err_bi_vec_pm30(i_8)', err_bi_vec_pm30(i_9)', err_bi_vec_pm50(i_35p)', err_bi_vec_pm50(i_40p)', err_bi_vec_pm50(i_45p)', err_bi_vec_pm50(i_50p)'];

err_low_plot_pm50 = [err_low_vec_pm50(i_50n)', err_low_vec_pm50(i_45n)', err_low_vec_pm50(i_40n)', ...
    err_low_vec_pm50(i_35n)', err_low_vec_pm30(i_1)', err_low_vec_pm30(i_2)', err_low_vec_pm15_20(i_1)', err_low_vec_pm15_20(i_2)', err_low_vec_pm10(i_3)',...
    err_low_vec_pm5(i_4)', ones(length(i_4),1)*low_nom, err_low_vec_pm5(i_6)', ...
    err_low_vec_pm10(i_7)', err_bi_vec_pm15_20(i_8)', err_bi_vec_pm15_20(i_9)',...
    err_low_vec_pm30(i_8)', err_low_vec_pm30(i_9)', err_low_vec_pm50(i_35p)', err_low_vec_pm50(i_40p)', err_low_vec_pm50(i_45p)', err_low_vec_pm50(i_50p)'];



delta_vec_pm50 = -50:5:50;

% save('results_plot/plot_data_pm50_r6', 'bound_plot_pm50', 'bi_plot_pm50', 'delta_vec_pm50', 'err_low_plot_pm50');

figure
subplot(1,2,1)
p1 = plot(delta_vec_pm50, 100*bound_plot_pm50, 'o-','LineWidth',LW,'MarkerSize',MS); 
xlabel('Delta $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); box on

subplot(1,2,2)
p1 = plot(delta_vec_pm50, 100*bi_plot_pm50, 'o-','LineWidth',LW,'MarkerSize',MS); 
xlabel('Delta $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('BF Error $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1(1),p1(2), p1(3), p1(4), p1(5), p1(6)],{'$\rho$','m', 'p', 't', 'a','SA1'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_2)

if save_on ==1
    saveas(gcf,'plots_airfoil/airfoil_line_pm50_r10','epsc')
    saveas(gcf,'plots_airfoil/airfoil_line_pm50_r10','png')
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 1000 samples
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% design_params_1000 = importdata(strcat('batch_1000_all', '/design_params.dat'), ',' , 1); 
% design_params_1000 = design_params_1000.data; 
% 
% load('batch_1000_res', 'cp_results_1000', 'err_bound_vec_1000', 'err_bi_vec_1000','err_low_vec_1000', 'N', 'n', 'r'); 
% 
% [bound_order, bound_i] = sort(err_bound_vec_1000); 
% [bi_order, bi_i] = sort(err_bi_vec_1000); 
% 
% figure
% p1 = plot(ones(1000,1)*100*bound_nom, '-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% hold on
% p2 = plot(ones(1000,1)*100*bi_nom, '--','Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(100*err_bound_vec_1000(bound_i), 'o','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
% p4 = plot(100*err_bi_vec_1000(bound_i), 'x','Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% legend([p3,p4, p1, p2],{'Bound Samples','BF Samples', 'Nominal Bound', 'Nominal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
% hold off
% xlabel('Sample $i$','interpreter','latex','Fontsize',FS)
% ylabel('Error $[\%]$','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% 
% if save_on ==1
%     saveas(gcf,'plots_airfoil/airfoil_rand_samples_all_r10','epsc')
%     saveas(gcf,'plots_airfoil/airfoil_rand_samples_all_r10','png')
% end
% 
% %%% Closer look - best 100 samples
% figure
% p1 = plot(100*ones(100,1)*bound_nom, '-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% hold on
% p2 = plot(100*ones(100,1)*bi_nom, '--','Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(100*err_bound_vec_1000(bound_i(1:100)), 'o','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
% p4 = plot(100*err_bi_vec_1000(bound_i(1:100)), 'x','Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% legend([p3,p4, p1, p2],{'Bound Samples','BF Samples', 'Nominal Bound', 'Nominal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
% hold off
% xlabel('Sample $i$','interpreter','latex','Fontsize',FS)
% ylabel('Error $[\%]$','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% 
% if save_on ==1
%     saveas(gcf,'plots_airfoil/airfoil_rand_samples_100_r10','epsc')
%     saveas(gcf,'plots_airfoil/airfoil_rand_samples_100_r10','png')
% end
% 
% % Observations:
% % Noisy data is not great. 
% 
% % % Check that bound is less then the best of +- 50. It is. 
% % err_bound_vec_1000(bound_i(1))
% % % vs 
% % bound_plot_pm50
% % and the difference in bi
% [~, index_min] = min(bound_plot_pm50(:));
% bi_temp = bi_plot_pm50(:); 
% bi_temp(index_min)
% 
% 
% %%% Table 
% % Best performing bound
% res_mat = 100*[low_nom, bound_nom, bi_nom; ...
%     err_low_vec_1000(bound_i(1)), err_bound_vec_1000(bound_i(1)), err_bi_vec_1000(bound_i(1)); ...
%     err_low_vec_1000(bi_i(1)), err_bound_vec_1000(bi_i(1)), err_bi_vec_1000(bi_i(1))];
% 
% res_tab = array2table(res_mat, 'VariableNames', {'LF', 'Bound', 'BF'}, 'RowNames', {'Nominal', 'Opt Bound', 'Opt BF'}); 
% res_tab
% 
% 
% % Best performing bi is much better than best performing bound.
% % Reduction from 7.6 to 4.5 or 1.3... % Want the 1.3
% % How to get the 1.3? Or closer to it? 
% % Look at best performing bound metrix... 
% % Challenging because best 3 spots are all not great. So noisy... Hmm
% 
% % Could say: 
% % 1. It doesn't work so well
% % 2. The correlation between bound and bi is noisy. Can contrast to the
% % L-shape
% 
% % Or 
% % 1. Inspect sample points where bound performed well/poor. Sample
% % different limits? This feels like cheating. 
% 
% % Create plotting script. 
% % Plot both types of best... Consider how to advance. I should probably
% % pause this and polish basis reduction... 
% 
% % Look at samples of best performing: 
% 
% design_params_1000(bound_i(1),1:6)
% design_params_1000(bi_i(1),1:6)
% 
% % density, m, p, t, a, sa1 - 7
% % Doesn't seem to be much to see here... 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% QoI - best
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Plot QoI for indices
% 
% cp_opt_bound_l = cp_results_1000(:,:,bound_i(1)); 
% cp_opt_bi_l = cp_results_1000(:,:,bi_i(1)); 
% cp_bound_nom = cp_results_pm5(:,:,1); 
% 
% %%%
% % BF estimate 
% rand_sample = 1:n; 
% 
% B_bound = cp_opt_bound_l/norm(cp_opt_bound_l,'fro');
% A = Uf/norm(Uf,'fro');
% [P_s,ix] = matrixIDvR(B_bound,r);
% cp_bi_bound = A(:,ix)*P_s*norm(Uf,'fro'); 
% 
% B_bi = cp_opt_bi_l/norm(cp_opt_bi_l,'fro');
% [P_s,ix] = matrixIDvR(B_bi,r);
% cp_bi_bi = A(:,ix)*P_s*norm(Uf,'fro'); 
% 
% B_bi = cp_bound_nom/norm(cp_bound_nom,'fro');
% [P_s,ix] = matrixIDvR(B_bi,r);
% cp_bi_nom = A(:,ix)*P_s*norm(Uf,'fro'); 
% 
% % % check bi estimate
% % norm(Uf-cp_bi_nom)/norm(Uf)
% % norm(Uf-cp_bi_bound)/norm(Uf)
% % norm(Uf-cp_bi_bi)/norm(Uf)
% 
% % %%%
% % Uf 
% % Uc_nom
% 
% % Locate run with largest imporovement in error. 
% [~, index_bound] = min(vecnorm(cp_bi_bound - Uf) - vecnorm(cp_bi_nom - Uf)); 
% [~, index_bi] = min(vecnorm(cp_bi_bi - Uf)- vecnorm(cp_bi_nom - Uf));
% 
% %%% Best bound 
% n_points = 200; 
% xx = linspace(-1,1,n_points);
% 
% figure
% p1 = plot(xx,Uf(:,index_bound),'-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% hold on
% p2 = plot(xx,Uc_nom(:,index_bound),'-', 'Color',c2,'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(xx,cp_bi_nom(:,index_bound),'--', 'Color',c3,'LineWidth',LW,'MarkerSize',MS); 
% p4 = plot(xx,cp_bi_bound(:,index_bound),'.-', 'Color',c4,'LineWidth',LW,'MarkerSize',MS); 
% hold off
% xlabel('$x$','interpreter','latex','Fontsize',FS)
% ylabel('Cp','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% legend([p1,p2, p3, p4],{'h','l', 'Nom', 'Opt'},'interpreter', 'latex', 'fontsize', FS_leg)
% 
% if save_on ==1
%     saveas(gcf,'plots_airfoil/airfoil_cp_realization_opt_bound_r10','epsc')
%     saveas(gcf,'plots_airfoil/airfoil_cp_realization_opt_bound_r10','png')
% end
% 
% figure
% p1 = plot(xx,Uf(:,index_bi),'-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% hold on
% p2 = plot(xx,Uc_nom(:,index_bi),'-', 'Color',c2,'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(xx,cp_bi_nom(:,index_bi),'--', 'Color',c3,'LineWidth',LW,'MarkerSize',MS); 
% p4 = plot(xx,cp_bi_bi(:,index_bi),':', 'Color',c4,'LineWidth',LW,'MarkerSize',MS); 
% hold off
% xlabel('$x$','interpreter','latex','Fontsize',FS)
% ylabel('Cp','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% legend([p1,p2, p3, p4],{'h','l', 'Nom', 'Opt'},'interpreter', 'latex', 'fontsize', FS_leg)
% 
% if save_on ==1
%     saveas(gcf,'plots_airfoil/airfoil_cp_realization_opt_bi_r10','epsc')
%     saveas(gcf,'plots_airfoil/airfoil_cp_realization_opt_bi_r10','png')
% end
% 
% 
% figure
% p1 = plot(xx,Uf(:,index_bound),'-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% hold on
% p2 = plot(xx,Uc_nom(:,index_bound),'-', 'Color',c2,'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(xx,cp_bi_nom(:,index_bound),'--', 'Color',c3,'LineWidth',LW,'MarkerSize',MS); 
% p4 = plot(xx,cp_bi_bound(:,index_bound),':', 'Color',c4,'LineWidth',LW,'MarkerSize',MS); 
% hold off
% xlabel('$x$','interpreter','latex','Fontsize',FS)
% ylabel('Cp','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% legend([p1,p2, p3, p4],{'h','l', 'Nom', 'Opt'},'interpreter', 'latex', 'fontsize', FS_leg)
% 
% if save_on ==1
%     saveas(gcf,'plots_airfoil/airfoil_cp_realization_opt_bound_r10','epsc')
%     saveas(gcf,'plots_airfoil/airfoil_cp_realization_opt_bound_r10','png')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% QoI - 10 best and 10 worst? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Histogram 1000
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Vizualize  data
% n_hist = 20; 
% 
% error_b_nom = vecnorm(cp_bi_nom - Uf)./vecnorm(Uf);
% error_b_opt = vecnorm(cp_bi_bound - Uf)./vecnorm(Uf);
% error_b_bi = vecnorm(cp_bi_bi - Uf)./vecnorm(Uf);
% 
% % check data on this one... 
% figure
% hold on
% l1 = xline(res_mat(1,2),'LineWidth', LW, 'Color', c1);
% l2 = xline(res_mat(2,2),'--','LineWidth', LW, 'Color', c2);
% 
% h1 = histogram(abs(100*error_b_nom),n_hist,'FaceColor',c1);
% h2 = histogram(abs(100*error_b_opt),n_hist,'FaceColor',c2);
% hold off
% legend([l1, l2, h1,h2],{'Nominal Bound','Optimal Bound','Nominal Ensemble','Optimal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
% xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
% ylabel('Frequency','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% set(gcf,'Position',size_1)
% 
% if save_on ==1
%     saveas(gcf,'plots_airfoil/airfoil_cp_hist_opt_bound_r10','epsc')
%     saveas(gcf,'plots_airfoil/airfoil_cp_hist_opt_bound_r10','png')
% end
% 
% figure
% hold on
% l1 = xline(res_mat(1,2),'LineWidth', LW, 'Color', c1);
% l2 = xline(res_mat(2,2),'--','LineWidth', LW, 'Color', c2);
% 
% h1 = histogram(abs(100*error_b_nom),n_hist,'FaceColor',c1);
% h2 = histogram(abs(100*error_b_bi),n_hist,'FaceColor',c2);
% hold off
% legend([l1, l2, h1,h2],{'Nominal Bound','Optimal Bound','Nominal Ensemble','Optimal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
% xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
% ylabel('Frequency','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% set(gcf,'Position',size_1)
% 
% if save_on ==1
%     saveas(gcf,'plots_airfoil/airfoil_cp_hist_opt_bi_r10','epsc')
%     saveas(gcf,'plots_airfoil/airfoil_cp_hist_opt_bi_r10','png')
% end
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Airfoil  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Plot optimal (with changed SA and density too?) 
% 
% %%% Load nominal coordinates
% folder_name_loop = 'limit_tests_pm5'; 
% i_ens = 1; 
% load(strcat(folder_name_loop, '/airfoil_x_data_', num2str(i_ens),'.mat'))
% load(strcat(folder_name_loop, '/airfoil_y_data_', num2str(i_ens),'.mat'))
% 
% % Ensemble of 500 :/ 
% plot(airfoil_x_data, airfoil_y_data, 'x')
% % Not very useful. 
% % airfoil_x_data, airfoil_x_data
% 
% design_params_1000(bound_i(1),1:6)
% design_params_1000(bi_i(1),1:6)
% % I wonder if the sample space is 'bumpier' near the chosen sample? A way
% % to improve heuristic/suggest further development? Ie look at best 10 and
% % plot. 
% 
% % figure
% % p1 = plot(ones(1000,1)*100*bound_nom, '-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% % hold on
% % p2 = plot(ones(1000,1)*100*bi_nom, '--','Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% % p3 = plot(100*err_bound_vec_1000(bound_i), 'o','Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
% % p4 = plot(100*err_bi_vec_1000(bound_i), 'x','Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% % legend([p3,p4, p1, p2],{'Bound Samples','BF Samples', 'Nominal Bound', 'Nominal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
% % hold off
% % xlabel('Sample $i$','interpreter','latex','Fontsize',FS)
% % ylabel('Error $[\%]$','interpreter','latex','Fontsize',FS)
% % axis tight
% % set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); box on
% 
% 
% err_bound_vec_1000(bound_i(1:10))
% design_params_1000(bound_i(1:10),1:6)
% 
% figure
% plot(design_params_1000(bound_i(1:10),1:6)','x','LineWidth',LW,'MarkerSize',MS)
% hold on
% plot(median(design_params_1000(bound_i(1:10),1:6)), 'o','LineWidth',LW,'MarkerSize',MS)
% 
% % Can I take the most common value, ie median? Give the sample a score of
% % center? 
% % This is definitely bias - increasing number of samples in window you would approach the center of
% % sampling space... 
% 
% % Measure distance from median... 
% est_middle = median(design_params_1000(bound_i(1:10),1:6));
% mid_metric = vecnorm(design_params_1000(bound_i(1:10)) - est_middle'); 
% 
% mid_mat = [mid_metric', 100*err_bound_vec_1000(bound_i(1:10))', 100*err_bi_vec_1000(bound_i(1:10))'];
% 
% [~, i_mid_sort] = sort(mid_mat(:,1)); 
% % If I order by middle? 
% mid_mat_sort = mid_mat(i_mid_sort, :); 
% 
% % Well that failed... 
% % Is there a different method to make a more conservative estimate from the
% % best handful of samples? 
% % It would make sense to use the info that the bound is noisy as a reason
% % not to trust just one point... unlike in the L-shaped elasticity.
% 
% % No! - this wouldn't make sense at all, because in practice this noise is
% % not know. This is a problem. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 500 samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


design_params_500 = importdata(strcat('batch_500_airfoil', '/design_params.dat'), ',' , 1); 
design_params_500 = design_params_500.data; 

% load('batch_500_res_r10', 'cp_results_500', 'err_bound_vec_500', 'err_bi_vec_500','err_low_vec_500', 'N', 'n', 'r'); 
load('batch_500_res', 'cp_results_500', 'err_bound_vec_500', 'err_bi_vec_500','err_low_vec_500', 'N', 'n', 'r'); 

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
    saveas(gcf,'plots_airfoil/airfoil_rand_samples_500_r6','epsc')
    saveas(gcf,'plots_airfoil/airfoil_rand_samples_500_r6','png')
end

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

if save_on ==1
    saveas(gcf,'plots_airfoil/airfoil_rand_samples_500_100_r6','epsc')
    saveas(gcf,'plots_airfoil/airfoil_rand_samples_500_100_r6','png')
end

res_mat_500 = 100*[low_nom, bound_nom, bi_nom; ...
    err_low_vec_500(bound_i(1)), err_bound_vec_500(bound_i(1)), err_bi_vec_500(bound_i(1)); ...
    err_low_vec_500(bi_i(1)), err_bound_vec_500(bi_i(1)), err_bi_vec_500(bi_i(1))];

res_tab_500 = array2table(res_mat_500, 'VariableNames', {'LF', 'Bound', 'BF'}, 'RowNames', {'Nominal', 'Opt Bound', 'Opt BF'}); 
res_tab_500

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Histogram 500
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cp_opt_bound_l = cp_results_500(:,:,bound_i(1)); 
cp_opt_bi_l = cp_results_500(:,:,bi_i(1)); 
cp_bound_nom = cp_results_pm5(:,:,1); 

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
    saveas(gcf,'plots_paper/airfoil_cp_hist_opt_bound_500_r6','epsc')
    saveas(gcf,'plots_airfoil/airfoil_cp_hist_opt_bound_500_r6','png')
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
    saveas(gcf,'plots_airfoil/airfoil_cp_hist_opt_bi_500_r6','epsc')
    saveas(gcf,'plots_airfoil/airfoil_cp_hist_opt_bi_500_r6','png')
end
