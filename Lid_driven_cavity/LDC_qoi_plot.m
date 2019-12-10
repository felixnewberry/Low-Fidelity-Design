

clear all
close all
clc

save_on = 0; 

% save_folder = 'plots_qoi';
% save_folder = 'plots_qoi_nump5';
% save_folder = 'plots_qoi_nums';
save_folder = 'plots_qoi_nups';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;     % Line width
MS = 8;     % Marker Size
FS_leg = 16; % Font size legend


size_1 = [0,0,445,345]; 
size_2 = [0,0,1340,515]; 

size_square = [0,0,445,445]; 
size_large = [0,0,668,518]; 

FS = 28;    % Font size axis
FS_axis = 18; 
LW_axis = 1; 

MyMarkerSize = 8; 

% Colors
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; 
c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880]; 
c6 = [0.3010, 0.7450, 0.9330]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% U Mid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('u_meshes/u_64_f_2.mat')

% 200 samples. 65 points
Uf_u = u_matrix_0'; 


% best u mid is at u 1.5579, 0.2737
% best p base (not found with routine bound method) is u +0.2175 nu +1.8316

% Plot velocity through mid plane 
% Maybe do for completeness

% Plot pressure along top

% need to change nominal and optimal rv 2. Then plot this again. 
% Difference will still likely be difficult to see. 

% plot the pressure differce too just to find out? Check which pressure. 

Uc_nom_u = load('LDC_design/u_mid_nom_2.mat', 'Uc','Ub','sb');
Ub_nom_u = Uc_nom_u.Ub; 
sb_nom_u = Uc_nom_u.sb; 
Uc_nom_u = Uc_nom_u.Uc; 

Uc_opt_u = load('LDC_design/u_mid_opt_2.mat', 'Uc','Ub','sb');
Ub_opt_u = Uc_opt_u.Ub; 
sb_opt_u = Uc_opt_u.sb; 
Uc_opt_u = Uc_opt_u.Uc; 

error_b_nom_u = vecnorm(Ub_nom_u-Uf_u)./vecnorm(Uf_u);
error_b_opt_u = vecnorm(Ub_opt_u-Uf_u)./vecnorm(Uf_u);

[~, index_max_u] = max(abs(error_b_nom_u - error_b_opt_u)); 
norm(Ub_opt_u(:,index_max_u) - Uf_u(:,index_max_u))
norm(Ub_nom_u(:,index_max_u) - Uf_u(:,index_max_u))

1; 

load 'x_64.mat'

x_highfidelity = x_64(:,1); 

plot_index = index_max_u; 

% % plot H, Nominal L, Nominal B
% plot u or v? 
% plot_uv = 1:33; 
% plot_uv = 34:66; 




figure 
p1 = plot(x_highfidelity,Uf_u(:,plot_index),'color',c1,'LineWidth',LW);
hold on
p2 = plot(x_highfidelity,Uc_nom_u(:,plot_index),':','color',c2,'LineWidth',LW);
p3 = plot(x_highfidelity,Ub_nom_u(:,plot_index),'-.','color',c3,'LineWidth',LW);
p4 = plot(x_highfidelity,Uc_opt_u(:,plot_index),':','color',c4,'LineWidth',LW);
p5 = plot(x_highfidelity,Ub_opt_u(:,plot_index),'-.','color',c5,'LineWidth',LW);
hold off
xlabel('$x$','interpreter','latex','Fontsize',FS)
ylabel('$ U_{mid}$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1,p2,p3,p4,p5],{'H','Nominal L','Nominal B', 'Optimal L', 'Optimal B'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3,p4,p5],{'H','Nominal L','Nominal B', 'Optimal L', 'Optimal B'},'interpreter', 'latex', 'fontsize', FS_leg)
title('Mid Velocity Realization','Interpreter','latex')
set(gcf,'Position',size_large)

if save_on ==1
    saveas(gcf,'plots_qoi/LDC_U_mid_realization','png')
end
1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% P Base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot base pressure

Uf_pb = u_matrix_3'; 

Uc_nom_pb = load('LDC_design/P_base_nom_2.mat', 'Uc','Ub','sb');
Ub_nom_pb = Uc_nom_pb.Ub; 
sb_nom_pb = Uc_nom_pb.sb; 
Uc_nom_pb = Uc_nom_pb.Uc; 

Uc_opt_pb = load('LDC_design/P_base_opt_2.mat', 'Uc','Ub','sb');
Ub_opt_pb = Uc_opt_pb.Ub; 
sb_opt_pb = Uc_opt_pb.sb; 
Uc_opt_pb = Uc_opt_pb.Uc; 

error_b_nom_pb = vecnorm(Ub_nom_pb-Uf_pb)./vecnorm(Uf_pb);
error_b_opt_pb = vecnorm(Ub_opt_pb-Uf_pb)./vecnorm(Uf_pb);

[~, index_max_pb] = max(abs(error_b_nom_pb - error_b_opt_pb)); 
norm(Ub_opt_pb(:,index_max_pb) - Uf_pb(:,index_max_pb))
norm(Ub_nom_pb(:,index_max_pb) - Uf_pb(:,index_max_pb))

load 'x_64.mat'

x_highfidelity = x_64(:,1); 

plot_index = index_max_pb; 

1; 




figure 
p1 = plot(x_highfidelity,Uf_pb(:,plot_index),'color',c1,'LineWidth',LW);
hold on
p2 = plot(x_highfidelity,Uc_nom_pb(:,plot_index),':','color',c2,'LineWidth',LW);
p3 = plot(x_highfidelity,Ub_nom_pb(:,plot_index),'-.','color',c3,'LineWidth',LW);
p4 = plot(x_highfidelity,Uc_opt_pb(:,plot_index),':','color',c4,'LineWidth',LW);
p5 = plot(x_highfidelity,Ub_opt_pb(:,plot_index),'-.','color',c5,'LineWidth',LW);
hold off
xlabel('$x$','interpreter','latex','Fontsize',FS)
ylabel('$P_{base} $','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1,p2,p3,p4,p5],{'H','Nominal L','Nominal B', 'Optimal L', 'Optimal B'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3,p4,p5],{'H','Nominal L','Nominal B', 'Optimal L', 'Optimal B'},'interpreter', 'latex', 'fontsize', FS_leg)
title('Base Pressure Realization','Interpreter','latex')
set(gcf,'Position',size_large)

if save_on ==1
    saveas(gcf,'plots_qoi/LDC_P_base_realization','png')
end

% first plot high-fidelity - then figure out the rest. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% U mag Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% center of vortex - for one sample... 

vortex_h = [0.6215, 0.7357]; 
vortex_l = [0.6215, 0.7357]; 
% they agree

y_lim_index = 65; 

step_quiv = 2; 

load('u_meshes/u_matrix_high.mat')

p_field_h = u_matrix_5; 
u_field_x_h = u_matrix_6; 
u_field_y_h = u_matrix_7; 

% load('u_meshes/u_matrix_low.mat')
load('u_meshes/u_matrix_low_nom.mat')


p_field_l = u_matrix_5; 
u_field_x_l = u_matrix_6; 
u_field_y_l = u_matrix_7; 


% Compare errors of mean or one individual run?

% I should first calculate errors and then the mean. 
p_mean_h = mean(p_field_h);
u_mean_x_h = mean(u_field_x_h);
u_mean_y_h = mean(u_field_y_h);

p_mean_l = mean(p_field_l);
u_mean_x_l = mean(u_field_x_l);
u_mean_y_l = mean(u_field_y_l);

nx = 64; 
% for now do 1st sample. 

x = x_64(:,1); 
y = x; 
% should I use contour - filled? yeah. 

% plot mean 
p_mat_h = reshape(p_mean_h,nx+1,nx+1); 
u_mat_h = reshape(u_mean_x_h,nx+1,nx+1); 
v_mat_h = reshape(u_mean_y_h,nx+1,nx+1); 
u_mag_h = sqrt(u_mat_h.^2 +v_mat_h.^2);

p_mat_l = reshape(p_mean_l,nx+1,nx+1); 
u_mat_l = reshape(u_mean_x_l,nx+1,nx+1); 
v_mat_l = reshape(u_mean_y_l,nx+1,nx+1); 
u_mag_l = sqrt(u_mat_l.^2 +v_mat_l.^2);

% u_min = min([u_mag_h(:);u_mag_l(:)]);
% u_max = max([u_mag_h(:);u_mag_l(:)]);
% 
% p_min = min([p_mat_h(:);p_mat_l(:)]);
% p_max = max([p_mat_h(:);p_mat_l(:)]);

% Plot mean velocity magnitude 


% 'LevelList',[-10,linspace(-0.1031, 0.0800,10),10]
% high fidelity

figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index,:),u_mag_l(1:y_lim_index,:),10);
set(h, 'edgecolor','none');
hold on
caxis_vec = caxis;
% caxis([u_min, u_max])
quiver(x(1:step_quiv:end),y(1:step_quiv:end),u_mat_l(1:step_quiv:end,1:...
    step_quiv:end), v_mat_l(1:step_quiv:end,1:step_quiv:end),'Color',c3)
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
pbaspect([1 1 1])
colorbar 
title('$U_L$','Interpreter', 'latex')


% low fidelity 
subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index,:),u_mag_h(1:y_lim_index,:));
set(h, 'edgecolor','none');
hold on
caxis(caxis_vec)
quiver(x(1:step_quiv:end),y(1:step_quiv:end),u_mat_h(1:step_quiv:end,1:...
    step_quiv:end), v_mat_h(1:step_quiv:end,1:step_quiv:end),'Color',c3)
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_2)
pbaspect([1 1 1])
colorbar 
title('$U_H$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,strcat(save_folder,'/LDC_U_mag'),'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% u field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

caxis_vec = [-0.1979    0.2403]; 
level_vec = [100*caxis_vec(1),linspace(caxis_vec(1), caxis_vec(2),9),caxis_vec(2)*100];

figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index,:),u_mat_l(1:y_lim_index,:),'LevelList',level_vec);
set(h, 'edgecolor','none');
hold on
caxis(caxis_vec)
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
pbaspect([1 1 1])
colorbar 
title('$u_L$','Interpreter', 'latex')


% low fidelity 
subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index,:),u_mat_h(1:y_lim_index,:),'LevelList',level_vec);
set(h, 'edgecolor','none');
hold on
caxis(caxis_vec)
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_2)
pbaspect([1 1 1])
colorbar 
title('$u_H$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,strcat(save_folder,'/LDC_u_field'),'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% v field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

caxis_vec = [-0.3534,   0.1359]; 
level_vec = [100*caxis_vec(1),linspace(caxis_vec(1), caxis_vec(2),9),caxis_vec(2)*100];

figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index,:),v_mat_l(1:y_lim_index,:),'LevelList',level_vec);
set(h, 'edgecolor','none');
hold on
caxis(caxis_vec)
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
pbaspect([1 1 1])
colorbar 
title('$v_L$','Interpreter', 'latex')

% low fidelity 
subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index,:),v_mat_h(1:y_lim_index,:),'LevelList',level_vec);
set(h, 'edgecolor','none');
hold on
caxis(caxis_vec)
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_2)
pbaspect([1 1 1])
colorbar 
title('$v_H$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,strcat(save_folder,'/LDC_v_field'),'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% P Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
caxis_vec = [-0.1, 0.1];
level_vec = [100*caxis_vec(1),linspace(caxis_vec(1), caxis_vec(2),9),caxis_vec(2)*100];

figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index,:),p_mat_l(1:y_lim_index,:),'LevelList',level_vec);
set(h, 'edgecolor','none');
caxis(caxis_vec)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
pbaspect([1 1 1])
colorbar
title('$P_L$','Interpreter', 'latex')


subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index,:),p_mat_h(1:y_lim_index,:),'LevelList',level_vec);
set(h, 'edgecolor','none');
caxis(caxis_vec)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_2)
pbaspect([1 1 1])
colorbar
title('$P_H$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,strcat(save_folder,'/LDC_P'),'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error - u, v, p and vort 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do l - h so positive means l overestimates. 

% l - h absolute error
error_u = u_field_x_l - u_field_x_h; 
error_v = u_field_y_l - u_field_y_h; 
error_p = p_field_l - p_field_h; 

vort_h = curl(u_field_x_h,u_field_y_h);
vort_l = curl(u_field_x_l,u_field_y_l);
error_vort = vort_l- vort_h; 

% mean and sd of 200 samples
e_mat_u_mean = reshape(mean(error_u),nx+1,nx+1); 
e_mat_v_mean = reshape(mean(error_v),nx+1,nx+1); 
e_mat_p_mean = reshape(mean(error_p),nx+1,nx+1); 

e_mat_u_sd = reshape(std(error_u),nx+1,nx+1); 
e_mat_v_sd = reshape(std(error_v),nx+1,nx+1); 
e_mat_p_sd = reshape(std(error_p),nx+1,nx+1); 

% Plot mean vorticity

% % compute curl and angular velocity from mean velocity
% [CURLZ_h, CAV_h] = curl(u_mat_h,v_mat_h);
% [CURLZ_l, CAV_l] = curl(u_mat_l,v_mat_l);

e_mat_vort_mean = reshape(mean(error_vort),nx+1,nx+1); 
e_mat_vort_sd = reshape(std(error_vort),nx+1,nx+1); 

vort_mean_h = reshape(mean(vort_h),nx+1,nx+1);
vort_mean_l = reshape(mean(vort_l),nx+1,nx+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% u velocity Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adjust limits to look at body of flow: 
% ie -0.1 to 0.1 and 0 to 0.03

caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

caxis_vec_sd = [0, 0.03];
level_vec_sd = linspace(caxis_vec_sd(1), caxis_vec_sd(2),11);

figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_square)
% pbaspect([1 1 1])
colorbar
title('u absolute error $\mu$','Interpreter', 'latex')


subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_sd(1:y_lim_index,:),'LevelList',level_vec_sd);
set(h, 'edgecolor','none');
caxis(caxis_vec_sd)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_square)
% pbaspect([1 1 1])
colorbar
set(gcf,'Position',size_2)
title('u absolute error $\sigma$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,strcat(save_folder,'/LDC_u_error'),'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% v velocity Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

caxis_vec_mean = [-0.2, 0.2];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];


figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('v absolute error $\mu$','Interpreter', 'latex')

subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_sd(1:y_lim_index,:));
set(h, 'edgecolor','none');
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_2)
% pbaspect([1 1 1])
colorbar
title('v absolute error $\sigma$','Interpreter', 'latex')

set(gcf,'Position',size_2)

if save_on ==1
    saveas(gcf,strcat(save_folder,'/LDC_v_error'),'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pressure Error 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

caxis_vec_sd = [0, 0.02];
level_vec_sd = linspace(caxis_vec_sd(1), caxis_vec_sd(2),11);

figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('P absolute error $\mu$','Interpreter', 'latex')

subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_sd(1:y_lim_index,:),'LevelList',level_vec_sd);
set(h, 'edgecolor','none');
caxis(caxis_vec_sd)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_2)
% pbaspect([1 1 1])
colorbar
title('P absolute error $\sigma$','Interpreter', 'latex')


if save_on ==1
    saveas(gcf,strcat(save_folder,'/LDC_P_error_trun'),'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Vorticity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_lim_index=65; 

caxis_vec = [-0.015, 0.015];
level_vec = [100*caxis_vec(1),linspace(caxis_vec(1), caxis_vec(2),10)];


% make these contour plots consistent. 
figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index,:),vort_mean_l(1:y_lim_index,:),'LevelList',level_vec);
set(h, 'edgecolor','none');
hold on
caxis(caxis_vec)
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
% xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_2)
pbaspect([1 1 1])
colorbar 
title('$\omega_L$ ','Interpreter', 'latex')


subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index,:),vort_mean_h(1:y_lim_index,:),'LevelList',level_vec); 
set(h, 'edgecolor','none');
hold on
caxis(caxis_vec)
% quiver(x(1:step_quiv:end),y(1:step_quiv:end),u_mat_h(1:step_quiv:end,1:...
%     step_quiv:end), v_mat_h(1:step_quiv:end,1:step_quiv:end),'Color',c3)
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
% xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
pbaspect([1 1 1])
colorbar 
title('$\omega_H$ ','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,strcat(save_folder,'/LDC_vort'),'png')
end

% figure
% [c,h]=contourf(x,y(1:y_lim_index,:),CURLZ_l(1:y_lim_index,:)-CURLZ_h(1:y_lim_index,:));
% set(h, 'edgecolor','none');
% hold on
% % caxis([-0.1031    0.0800])
% % quiver(x(1:step_quiv:end),y(1:step_quiv:end),u_mat_h(1:step_quiv:end,1:...
% %     step_quiv:end), v_mat_h(1:step_quiv:end,1:step_quiv:end),'Color',c3)
% plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
% xlabel('x','interpreter', 'latex', 'fontsize', FS)
% ylabel('y','interpreter', 'latex', 'fontsize', FS)
% axis tight
% % xlim([0,1]); ylim([0,1]);
% new_labels = linspace(0, 1, 3);
% set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
% colorbar 
% title('$\omega_H$ ','Interpreter', 'latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Vorticity Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_lim_index = 65; 

caxis_vec_mean = [-0.015, 0.015];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

caxis_vec_sd = [0, 11e-3];
level_vec_sd = linspace(caxis_vec_sd(1), caxis_vec_sd(2),11);


figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_vort_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$\omega$ absolute error $\mu$','Interpreter', 'latex')

subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_vort_sd(1:y_lim_index,:),'LevelList',level_vec_sd);
set(h, 'edgecolor','none');
caxis(caxis_vec_sd)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_2)
% pbaspect([1 1 1])
colorbar
title('$\omega$ absolute error $\sigma$','Interpreter', 'latex')

set(gcf,'Position',size_2)

if save_on ==1
    saveas(gcf,strcat(save_folder,'/LDC_vort_error'),'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot to illustrate general effects for u,v and nu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vortex_h = [0.6215, 0.7357]; 

p_center = [0.6684, 0.7357]; 

% x p is 0.7138, y is 0.7357. 
% works for v but not so well for u 
% or 0.7138 for both



load('u_meshes/u_matrix_low_nom.mat') % nominal 

p_field_l = u_matrix_5; 
u_field_x_l = u_matrix_6; 
u_field_y_l = u_matrix_7; 

% l - h absolute error
error_u = u_field_x_l - u_field_x_h; 
error_v = u_field_y_l - u_field_y_h; 
error_p = p_field_l - p_field_h; 

% mean and sd of 200 samples
e_mat_u_mean = reshape(mean(error_u),nx+1,nx+1); 
e_mat_v_mean = reshape(mean(error_v),nx+1,nx+1); 
e_mat_p_mean = reshape(mean(error_p),nx+1,nx+1); 

figure

subplot(3,3,1)
caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
plot(p_center(1),p_center(2),'x','Color',c2,'MarkerSize',MS,'LineWidth',LW)

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$u$ nominal $\nu$','Interpreter', 'latex')


caxis_vec_mean = [-0.2, 0.2];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];


subplot(3,3,2)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
plot(p_center(1),p_center(2),'x','Color',c2,'MarkerSize',MS,'LineWidth',LW)

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$v$ nominal $\nu$','Interpreter', 'latex')

subplot(3,3,3)
caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
plot(p_center(1),p_center(2),'x','Color',c2,'MarkerSize',MS,'LineWidth',LW)

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$P$ nominal $\nu$','Interpreter', 'latex')

% - 65 and + 100 % nu 
% Then try pointwise Gaussian


load('u_meshes/u_matrix_low_mp65.mat') % - 65 %  

p_field_l = u_matrix_5; 
u_field_x_l = u_matrix_6; 
u_field_y_l = u_matrix_7; 

% l - h absolute error
error_u = u_field_x_l - u_field_x_h; 
error_v = u_field_y_l - u_field_y_h; 
error_p = p_field_l - p_field_h; 

% mean and sd of 200 samples
e_mat_u_mean = reshape(mean(error_u),nx+1,nx+1); 
e_mat_v_mean = reshape(mean(error_v),nx+1,nx+1); 
e_mat_p_mean = reshape(mean(error_p),nx+1,nx+1); 


subplot(3,3,4)

caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$u$ $\nu -65 \%$','Interpreter', 'latex')


caxis_vec_mean = [-0.2, 0.2];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];


subplot(3,3,5)

[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$v$ $\nu -65 \%$','Interpreter', 'latex')

subplot(3,3,6)
caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$P$ $\nu-65 \%$','Interpreter', 'latex')


load('u_meshes/u_matrix_low_p1.mat') % + 100 %  

p_field_l = u_matrix_5; 
u_field_x_l = u_matrix_6; 
u_field_y_l = u_matrix_7; 

% l - h absolute error
error_u = u_field_x_l - u_field_x_h; 
error_v = u_field_y_l - u_field_y_h; 
error_p = p_field_l - p_field_h; 

% mean and sd of 200 samples
e_mat_u_mean = reshape(mean(error_u),nx+1,nx+1); 
e_mat_v_mean = reshape(mean(error_v),nx+1,nx+1); 
e_mat_p_mean = reshape(mean(error_p),nx+1,nx+1); 

subplot(3,3,7)

caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$u$ $\nu+100 \%$','Interpreter', 'latex')


caxis_vec_mean = [-0.2, 0.2];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];


subplot(3,3,8)

[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$v$ $\nu+100 \%$','Interpreter', 'latex')

subplot(3,3,9)
caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$P$ $\nu+100 \%$','Interpreter', 'latex')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot to illustrate sigmoid effects for u,v and nu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vortex_h = [0.6215, 0.7357]; 

p_center = [0.6684, 0.7357]; 

% x p is 0.7138, y is 0.7357. 
% works for v but not so well for u 
% or 0.7138 for both



load('u_meshes/u_matrix_low_nom.mat') % nominal 

p_field_l = u_matrix_5; 
u_field_x_l = u_matrix_6; 
u_field_y_l = u_matrix_7; 

% l - h absolute error
error_u = u_field_x_l - u_field_x_h; 
error_v = u_field_y_l - u_field_y_h; 
error_p = p_field_l - p_field_h; 

% mean and sd of 200 samples
e_mat_u_mean = reshape(mean(error_u),nx+1,nx+1); 
e_mat_v_mean = reshape(mean(error_v),nx+1,nx+1); 
e_mat_p_mean = reshape(mean(error_p),nx+1,nx+1); 

figure

subplot(3,3,1)
caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
plot(p_center(1),p_center(2),'x','Color',c2,'MarkerSize',MS,'LineWidth',LW)

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$u$ nominal $\nu$','Interpreter', 'latex')


caxis_vec_mean = [-0.2, 0.2];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];


subplot(3,3,2)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
plot(p_center(1),p_center(2),'x','Color',c2,'MarkerSize',MS,'LineWidth',LW)

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$v$ nominal $\nu$','Interpreter', 'latex')

subplot(3,3,3)
caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
plot(p_center(1),p_center(2),'x','Color',c2,'MarkerSize',MS,'LineWidth',LW)

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$P$ nominal $\nu$','Interpreter', 'latex')

% - 65 and + 100 % nu 
% Then try pointwise Gaussian


load('u_meshes/u_matrix_low_s1_mp5.mat') % +100% surrounding, -50% at vortex

p_field_l = u_matrix_5; 
u_field_x_l = u_matrix_6; 
u_field_y_l = u_matrix_7; 

% l - h absolute error
error_u = u_field_x_l - u_field_x_h; 
error_v = u_field_y_l - u_field_y_h; 
error_p = p_field_l - p_field_h; 

% mean and sd of 200 samples
e_mat_u_mean = reshape(mean(error_u),nx+1,nx+1); 
e_mat_v_mean = reshape(mean(error_v),nx+1,nx+1); 
e_mat_p_mean = reshape(mean(error_p),nx+1,nx+1); 


subplot(3,3,4)

caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$u$ Sigmoid $\nu$ $100$ to $-50$  $\%$','Interpreter', 'latex')


caxis_vec_mean = [-0.2, 0.2];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];


subplot(3,3,5)

[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$v$ Sigmoid $\nu$ $100$ to $-50$  $\%$','Interpreter', 'latex')

subplot(3,3,6)
caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$P$ Sigmoid $\nu$ $100$ to $-50$  $\%$','Interpreter', 'latex')


load('u_meshes/u_matrix_low_mp5_s1.mat') % + 100 %  

p_field_l = u_matrix_5; 
u_field_x_l = u_matrix_6; 
u_field_y_l = u_matrix_7; 

% l - h absolute error
error_u = u_field_x_l - u_field_x_h; 
error_v = u_field_y_l - u_field_y_h; 
error_p = p_field_l - p_field_h; 

% mean and sd of 200 samples
e_mat_u_mean = reshape(mean(error_u),nx+1,nx+1); 
e_mat_v_mean = reshape(mean(error_v),nx+1,nx+1); 
e_mat_p_mean = reshape(mean(error_p),nx+1,nx+1); 

subplot(3,3,7)

caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$u$ Sigmoid $\nu$ $-50$ to $100$  $\%$','Interpreter', 'latex')


caxis_vec_mean = [-0.2, 0.2];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];


subplot(3,3,8)

[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$v$ Sigmoid $\nu$ $-50$ to $100$  $\%$','Interpreter', 'latex')

subplot(3,3,9)
caxis_vec_mean = [-0.1, 0.1];
level_vec_mean = [100*caxis_vec_mean(1),linspace(caxis_vec_mean(1), caxis_vec_mean(2),10)];

[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_mean(1:y_lim_index,:),'LevelList',level_vec_mean);
set(h, 'edgecolor','none');
caxis(caxis_vec_mean)
hold on
plot(vortex_h(1),vortex_h(2),'o','Color',c2,'MarkerSize',MS,'LineWidth',LW)
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% pbaspect([1 1 1])
colorbar
title('$P$ Sigmoid $\nu$ $-50$ to $100$  $\%$','Interpreter', 'latex')