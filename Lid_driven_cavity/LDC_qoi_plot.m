

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
%%% QoI Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


y_lim_index = 50; 

step_quiv = 2; 

load('u_meshes/u_matrix_high.mat')

p_field_h = u_matrix_5; 
u_field_x_h = u_matrix_6; 
u_field_y_h = u_matrix_7; 

load('u_meshes/u_matrix_low.mat')

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

% u_sd_h = std(u_field_x_h); 
% 
% % look at the mean and sd. 
% n_hist = 30; 
% 
% % some values are close to zero. 
% figure
% plot(u_mean_x_h, 'x', 'color', c1); 
% hold on 
% plot(u_sd_h,'o','color',c2); 
% 
% figure 
% plot(u_sd_h./u_mean_x_h,'x', 'color', c1)
% figure
% hold on
% h1 = histogram(100*u_mean_h,n_hist,'FaceColor',c1);
% h2 = histogram(100*u_sd_h,n_hist,'FaceColor',c2);
% hold off
% legend([h1,h2],{'mean','sd'},'interpreter', 'latex', 'fontsize', FS_leg)
% xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
% ylabel('Frequency','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% % grid on
% set(gcf,'Position',size_1)

nx = 64; 
% for now do 1st sample. 

x = x_64(:,1); 
y = x; 
% should I use contour - filled? yeah. 

% plot mean 
p_mat_h = reshape(p_mean_h(1,:),nx+1,nx+1); 
u_mat_h = reshape(u_mean_x_h(1,:),nx+1,nx+1); 
v_mat_h = reshape(u_mean_y_h(1,:),nx+1,nx+1); 
u_mag_h = sqrt(u_mat_h.^2 +v_mat_h.^2);

p_mat_l = reshape(p_mean_l(1,:),nx+1,nx+1); 
u_mat_l = reshape(u_mean_x_l(1,:),nx+1,nx+1); 
v_mat_l = reshape(u_mean_y_l(1,:),nx+1,nx+1); 
u_mag_l = sqrt(u_mat_l.^2 +v_mat_l.^2);

u_min = min([u_mag_h(:);u_mag_l(:)]);
u_max = max([u_mag_h(:);u_mag_l(:)]);

p_min = min([p_mat_h(:);p_mat_l(:)]);
p_max = max([p_mat_h(:);p_mat_l(:)]);

% Plot a sample flow 
figure
[c,h]=contourf(x,y(1:y_lim_index,:),u_mag_h(1:y_lim_index,:));
set(h, 'edgecolor','none');
hold on
% caxis([u_min, u_max])
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
set(gcf, 'Position', size_square)
pbaspect([1 1 1])
colorbar 

if save_on ==1
    saveas(gcf,'plots_qoi/LDC_U_field_sample_h_trun','png')
end

figure
[c,h]=contourf(x,y(1:y_lim_index,:),u_mag_l(1:y_lim_index,:));
set(h, 'edgecolor','none');
hold on
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
set(gcf, 'Position', size_square)
pbaspect([1 1 1])
colorbar 
if save_on ==1
    saveas(gcf,'plots_qoi/LDC_U_field_sample_l_trun','png')
end

1; 

% save contour? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% P Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure field looks silly when comparing. Take base pressure instead?
% Already have that... 
figure
[c,h]=contourf(x,y(1:y_lim_index,:),p_mat_h(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
pbaspect([1 1 1])
colorbar
if save_on ==1
    saveas(gcf,'plots_qoi/LDC_P_field_sample_h','png')
end

figure
[c,h]=contourf(x,y(1:y_lim_index,:),p_mat_l(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
pbaspect([1 1 1])
colorbar
if save_on ==1
    saveas(gcf,'plots_qoi/LDC_P_field_sample_l','png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error - u, v and p sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do l - h so positive means l overestimates. 



error_u = u_field_x_l - u_field_x_h; 
error_v = u_field_y_l - u_field_y_h; 
error_p = p_field_l - p_field_h; 

e_mat_u = reshape(error_u(1,:),nx+1,nx+1); 
e_mat_v = reshape(error_v(1,:),nx+1,nx+1); 
e_mat_p = reshape(error_p(1,:),nx+1,nx+1); 




% an error sample. If I instead plot the mean error and the sd, for u v and
% p. 
figure
[c,h]=contourf(x,y(1:y_lim_index),e_mat_u(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
pbaspect([1 1 1])

figure
[c,h]=contourf(x,y(1:y_lim_index),e_mat_v(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
pbaspect([1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error - u, v and p mean and sd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_mat_u_mean = reshape(mean(error_u),nx+1,nx+1); 
e_mat_v_mean = reshape(mean(error_v),nx+1,nx+1); 
e_mat_p_mean = reshape(mean(error_p),nx+1,nx+1); 

e_mat_u_sd = reshape(std(error_u),nx+1,nx+1); 
e_mat_v_sd = reshape(std(error_v),nx+1,nx+1); 
e_mat_p_sd = reshape(std(error_p),nx+1,nx+1); 

% u 
figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_mean(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
title('u $\mu$','Interpreter', 'latex')


subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_u_sd(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
title('u $\sigma$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots_qoi/LDC_u_error_trun','png')
end

figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_mean(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
title('v $\mu$','Interpreter', 'latex')

subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_v_sd(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
title('v $\sigma$','Interpreter', 'latex')

set(gcf,'Position',size_2)

if save_on ==1
    saveas(gcf,'plots_qoi/LDC_v_error_trun','png')
end

figure
subplot(1,2,1)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_mean(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
title('P $\mu$','Interpreter', 'latex')

subplot(1,2,2)
[c,h]=contourf(x,y(1:y_lim_index),e_mat_p_sd(1:y_lim_index,:));
set(h, 'edgecolor','none');
% caxis([p_min, p_max])
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
title('P $\sigma$','Interpreter', 'latex')

set(gcf,'Position',size_2)

if save_on ==1
    saveas(gcf,'plots_qoi/LDC_P_error_trun','png')
end

