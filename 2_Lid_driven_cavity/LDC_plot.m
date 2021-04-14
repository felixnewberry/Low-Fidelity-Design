% Felix Newberry
% Date: 08-28-19

% Lid driven cavity
% Plot PC results of samples 

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
c7 = [0.6350, 0.0780, 0.1840]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plot ghia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
load 'x_32.mat'
load 'y_32.mat'

v_32 = load('u_y_array_32.mat');
v_32 = v_32.u_y_array; 
u_32 = load('u_x_array_32.mat','u_x_array');
u_32 = u_32.u_x_array; 

v_8 = load('u_y_array_8.mat');
v_8 = v_8.u_y_array; 
u_8 = load('u_x_array_8.mat');
u_8 = u_8.u_x_array; 

v_6 = load('u_y_array_6.mat');
v_6 = v_6.u_y_array; 
u_6 = load('u_x_array_6.mat');
u_6 = u_6.u_x_array; 

v_4 = load('u_y_array_4.mat');
v_4 = v_4.u_y_array; 
u_4 = load('u_x_array_4.mat');
u_4 = u_4.u_x_array; 

% Ghiea et al
% Results for u velocity along Vertical Line through Geometric Center of
% Cavityu_out
u_ghia_100 = [1, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, ...
    -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, ...
    -0.04775, -0.04192, -0.03717, 0];
y_ghia =[1, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, ...
    0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0];

v_ghia_100 = [0, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, ...
    -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, ...
    0.10890, 0.10091, 0.09233, 0];
x_ghia = [1, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, ...
    0.5, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0];

figure
subplot(2,1,1)
hold on
p2 = plot(u_32, y_32(:,2),'Color',c1,'LineWidth',LW);
p3 = plot(u_8, y_32(:,2),'Color',c2,'LineWidth',LW);
p4 = plot(u_6, y_32(:,2),'Color',c3,'LineWidth',LW);
p5 = plot(u_4, y_32(:,2),'Color',c4,'LineWidth',LW);
p1 = plot(u_ghia_100, y_ghia,'ro','LineWidth',LW);

ylabel('$u$ velocity', 'interpreter', 'latex', 'fontsize', FS)
xlabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
%legend('ghia et al', 'fenics')
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
axis tight;
hold off

subplot(2,1,2)
hold on
p2 = plot(x_32(:,1),v_32, 'Color',c1,'LineWidth',LW);
p3 = plot(x_32(:,1),v_8, 'Color',c2,'LineWidth',LW);
p4 = plot(x_32(:,1),v_6, 'Color',c3,'LineWidth',LW);
p5 = plot(x_32(:,1),v_4, 'Color',c4,'LineWidth',LW);
p1 = plot(x_ghia, v_ghia_100,'ro','LineWidth',LW);

xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$v$ velocity', 'interpreter', 'latex', 'fontsize', FS)
legend([p1,p2,p3,p4,p5],{'ghia et al','32','8','6','4'},...
       'interpreter', 'latex', 'fontsize', FS_leg/2)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
axis tight;
set(gcf,'Position',size_1)




%%% Plot contour and then streamlines


load('LDC_data/stream_coords')
x_stream = x; 
y_stream = y; 

load 'LDC_data/u_out_plot.mat'

n_cell_x = 128; %CHANGE TO 128.
nx = n_cell_x;
x = linspace(0,1,n_cell_x+1);


x = (x-0.5)*2;
x = 0.5*(cos(pi*(x-1)/2)+1);
y = x;

u = u_out(1:end/2);
v = u_out(end/2+1:end);


u_mat = reshape(u,nx+1,nx+1)'; 
v_mat = reshape(v,nx+1,nx+1)'; 

u_mag = sqrt(u_mat.^2 +v_mat.^2);

x_1 = linspace(0,1.0,11);
y_1 = 0.0*ones(1, length(x_1)); 
y_2 = 0.5*ones(1, length(x_1)); 


step_quiv = 10; 

figure
[c,h]=contourf(x,y,u_mag);
set(h, 'edgecolor','none');
hold on

quiver(x(1:step_quiv:end),y(1:step_quiv:end),u_mat(1:step_quiv:end,1:step_quiv:end), v_mat(1:step_quiv:end,1:step_quiv:end),'Color',c3)

% plot qoi
p1 = plot(x_1, y_2,'-.','color',c2,'LineWidth',LW+1);
p2 = plot(x_1, y_1,'--','color',c3,'LineWidth',LW+1);

X_arrow = [0.35 0.7];
Y_arrow = [0.95   0.95];
hh = annotation('arrow',X_arrow,Y_arrow,'Color','r');
set(hh, 'LineWidth', LW)

xlabel('$x$','interpreter', 'latex', 'fontsize', FS)
ylabel('$y$','interpreter', 'latex', 'fontsize', FS)
legend([p1,p2],{'$U_{\mathrm{mid}}$','$P_{\mathrm{base}}$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis+5, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf, 'Position', size_square)
pbaspect([1 1 1])

if save_on ==1
    saveas(gcf,'plots/LDC_Geom_contour','epsc')
end

% try quiver... 
% select step size: 



% Option1, standard. 
figure
h = streamslice(x,y, u_mat, v_mat, 1, 'cubic'); %, 0.9,0.1, [0.1, 8000])
set( h, 'Color', c1 )
hold on
set(h, 'LineWidth', LW) %LW/2 ?

% plot qoi
p1 = plot(x_1, y_2,'-.','color',c2,'LineWidth',LW+1);
p2 = plot(x_1, y_1,'--','color',c3,'LineWidth',LW+1);

X_arrow = [0.35 0.7];
Y_arrow = [0.95   0.95];
hh = annotation('arrow',X_arrow,Y_arrow,'Color','r');
set(hh, 'LineWidth', LW)

hold off
xlabel('$x$','interpreter', 'latex', 'fontsize', FS)
ylabel('$y$','interpreter', 'latex', 'fontsize', FS)
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis+5, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% set(gcf, 'Position', size_1)
set(gcf, 'Position', size_square)
legend([p1,p2],{'$U_{\mathrm{mid}}$','$P_{\mathrm{base}}$'},'interpreter', 'latex', 'fontsize', FS_leg)

pbaspect([1 1 1])
if save_on ==1
    saveas(gcf,'plots/LDC_Geom_stream','epsc')
end

1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_l = 4; 
n_h = 64; 

x_l  = linspace(0,1,n_l); 
x_h  = linspace(0,1,n_h); 


% concentrate points close to sides
x_l = (x_l - 0.5).*2; 
x_l = 0.5.*(cos(pi.*(x_l-1) / 2) + 1);

[xx_l, yy_l ] = meshgrid(x_l); 

x_h = (x_h - 0.5).*2; 
x_h = 0.5.*(cos(pi.*(x_h-1) / 2) + 1);

[xx_h, yy_h ] = meshgrid(x_h); 

% plot mesh.. 
% plot
figure
plot(xx_l, yy_l, 'k','LineWidth',LW/2)
hold on; 
plot(yy_l, xx_l, 'k','LineWidth',LW/2)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
% grid on; 
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis+5, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
axis tight;

set(gcf,'Position',size_square)

if save_on ==1
    saveas(gcf,'plots/mesh_low','epsc')
end

figure
plot(xx_h, yy_h, 'k','LineWidth',LW/4)
hold on; 
plot(yy_h, xx_h, 'k','LineWidth',LW/4)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
% grid on; 
set(gca,'Fontsize', FS_axis+5, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
axis tight;
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gcf,'Position',size_square)


if save_on ==1
    saveas(gcf,'plots/mesh_high','epsc')
end

1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Individual Sensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% r = 1, n = 3 ... 


load('LDC_design/individual_qoi_nu.mat')
% load('LDC_design/line_qoi_nu_2_n8.mat')
% load('LDC_design/line_qoi_nu_1.mat')

delta_nu_vec = delta_vec; 
error_bound_nu = error_bound_mat; 

load('LDC_design/individual_qoi_u.mat')
% load('LDC_design/line_qoi_u_1.mat')
% load('LDC_design/line_qoi_u_2_n8.mat')

delta_u_vec = delta_vec; 
error_bound_u = error_bound_mat; 

% plot nu change for different QoI, ie u_mid, u_field, p_field, p_top
figure
hold on
p1 = plot(100*delta_nu_vec,100*error_bound_nu(1,:),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_nu_vec,100*error_bound_nu(2,:),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(100*delta_nu_vec,100*error_bound_nu(3,:),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
% p4 = plot(100*delta_nu_vec,100*error_bound_nu(4,:),'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
p5 = plot(100*delta_nu_vec,100*error_bound_nu(5,:),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$\Delta \nu [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([p1,p5],{'$U_{\mathrm{mid}}$','$P_{\mathrm{base}}$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
xlim([-50,300]); 
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity $\nu$','Interpreter', 'latex')

if save_on ==1
saveas(gcf,'plots/LDC_nu','epsc')
%     saveas(gcf,'plots/LDC_1_nu','epsc')
%     saveas(gcf,'Plots/LDC_nu','epsc')
%     saveas(gcf,'plots/LDC_2_nu_8','epsc')
end


% plot u change for different QoI, ie u_mid, u_field, p_field, p_top
figure
hold on
p1 = plot(100*delta_u_vec,100*error_bound_u(1,:),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_u_vec,100*error_bound_u(2,:),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(100*delta_u_vec,100*error_bound_u(3,:),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
% p4 = plot(100*delta_u_vec,100*error_bound_u(4,:),'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
p5 = plot(100*delta_u_vec,100*error_bound_u(5,:),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$\Delta U [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4, p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([p1 p5],{'$U_{\mathrm{mid}}$','$P_{\mathrm{base}}$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
xlim([-70,200]); 
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity $U$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/LDC_u','epsc')
end

1; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Grid search - u and vanilla nu
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('LDC_design/grid_search.mat')


%%% Plot response surface 
% Step through each qoi 
qoi_vec = 1:5; 

plot_label = ["$U_{\mathrm{mid}}$","$U_{\mathrm{vert}}$", "$P$ Mid", "$P$ Vert", "$P_{\mathrm{base}}$"]; 
plot_save = ["u_Mid","U_Vert", "P_Mid", "P_Vert", "P_Base"]; 


%%%%%%%%%%%%%%%%%%%
[xrow,ycol] = meshgrid(...
            linspace(100*delta_nu_vec(1),100*delta_nu_vec(end),20),...
            linspace(100*delta_u_vec(1),100*delta_u_vec(end),20)...
          );
[xq,yq] = meshgrid(...
            linspace(100*delta_nu_vec(1),100*delta_nu_vec(end),2000),...
            linspace(100*delta_u_vec(1),100*delta_u_vec(end),2000)...
          );
      
% exclude 0s - these are nans. 
error_bound_mat(error_bound_mat==0) = nan;

error_bound_matq = zeros(length(qoi_vec),2000,2000); 
error_Bi_matq = zeros(length(qoi_vec),2000,2000); 

min_bound = zeros(length(qoi_vec),1); 
% min_bound_i
u_bound = zeros(length(qoi_vec),1); 
nu_bound = zeros(length(qoi_vec),1); 

min_bi = zeros(length(qoi_vec),1); 
% min_bi_i
u_bi = zeros(length(qoi_vec),1); 
nu_bi = zeros(length(qoi_vec),1); 

min_bi_opt = zeros(length(qoi_vec),1); 
%%% Interpolate 
for i_qoi = 1:length(qoi_vec)
    error_bound_matq(i_qoi,:,:) = interp2(xrow,ycol,100*reshape(error_bound_mat(i_qoi,:,:),length(delta_nu_vec), length(delta_u_vec))',xq,yq,'linear');
    error_Bi_matq(i_qoi,:,:) = interp2(xrow,ycol,100*reshape(error_Bi_mat(i_qoi,:,:),length(delta_nu_vec), length(delta_u_vec))',xq,yq,'linear');
    
    [min_bound(i_qoi), min_bound_i] = min(reshape(error_bound_matq(i_qoi,:,:),[],1)); 
    u_bound(i_qoi) = yq(min_bound_i); 
    nu_bound(i_qoi) = xq(min_bound_i); 
    
    [min_bi(i_qoi), min_bi_i] = min(reshape(error_Bi_matq(i_qoi,:,:),[],1)); 
    u_bi(i_qoi) = yq(min_bi_i); 
    nu_bi(i_qoi) = xq(min_bi_i); 
    
    bi_temp = reshape(error_Bi_matq(i_qoi,:,:),[],1); 
    min_bi_opt(i_qoi) = bi_temp(min_bound_i); 
end

max_z = max(error_bound_matq(:)); 

for i_qoi = 1:length(qoi_vec)

figure
hold on
surf(xq,yq,reshape(error_bound_matq(i_qoi,:,:),length(xq),length(yq)))
shading interp
view(0,90)
colorbar
% 5, 0.2
% p1 = plot3(0,0,max_z,'o','Color',c2,'MarkerSize',8,'linewidth',LW);
% p3 = plot3(nu_bi(i_qoi),u_bi(i_qoi),max_z,'d','Color',c7,'MarkerSize',8,'linewidth',LW);
% p2 = plot3(nu_bound(i_qoi),u_bound(i_qoi),max_z,'k+','MarkerSize',8,'linewidth',LW);
p1 = plot3(0,0,max_z,'ro','MarkerSize',8,'linewidth',LW);
p3 = plot3(nu_bi(i_qoi),u_bi(i_qoi),max_z,'rd','MarkerSize',8,'linewidth',LW);
p2 = plot3(nu_bound(i_qoi),u_bound(i_qoi),max_z,'r+','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Nominal','Optimal', 'True Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('$\Delta \nu $ [\%]','interpreter','latex','Fontsize',FS)
ylabel('$\Delta u$ [\%]','interpreter','latex','Fontsize',FS)
axis tight
caxis([min(reshape([error_bound_matq(i_qoi,:,:); error_Bi_matq(i_qoi,:,:)],1,[])) max(reshape([error_bound_matq(i_qoi,:,:); error_Bi_matq(i_qoi,:,:)],1,[]))])
c =colorbar;
c.TickLabelInterpreter = 'latex'; 
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,strcat('plots/LDC_', plot_save(i_qoi), '_bound'),'epsc')
end


figure
hold on
surf(xq,yq,reshape(error_Bi_matq(i_qoi,:,:),length(xq),length(yq)))
shading interp
view(0,90)
colorbar
% 5, 0.2
% p1 = plot3(0,0,max_z,'o','Color',c2,'MarkerSize',8,'linewidth',LW);
% p3 = plot3(nu_bi(i_qoi),u_bi(i_qoi),max_z,'d','Color',c7,'MarkerSize',8,'linewidth',LW);
% p2 = plot3(nu_bound(i_qoi),u_bound(i_qoi),max_z,'k+','MarkerSize',8,'linewidth',LW);
p1 = plot3(0,0,max_z,'ro','MarkerSize',8,'linewidth',LW);
p3 = plot3(nu_bi(i_qoi),u_bi(i_qoi),max_z,'rd','MarkerSize',8,'linewidth',LW);
p2 = plot3(nu_bound(i_qoi),u_bound(i_qoi),max_z,'r+','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Nominal','Optimal', 'True Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('$\Delta \nu $ [\%]','interpreter','latex','Fontsize',FS)
ylabel('$\Delta u$ [\%]','interpreter','latex','Fontsize',FS)
axis tight
caxis([min(reshape([error_bound_matq(i_qoi,:,:); error_Bi_matq(i_qoi,:,:)],1,[])) max(reshape([error_bound_matq(i_qoi,:,:); error_Bi_matq(i_qoi,:,:)],1,[]))])
c =colorbar;
c.TickLabelInterpreter = 'latex'; 
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title(strcat(plot_label(i_qoi),' BF Error'),'Interpreter','latex')



if save_on ==1
        saveas(gcf,strcat('plots/LDC_', plot_save(i_qoi), '_bi'),'epsc')
end

 1;
 

end

optimal_mat = [u_bound'; nu_bound']; 

optimal_tab = array2table(optimal_mat,...
    'VariableNames',{'U_Mid', 'U_Vert', 'P_Mid', 'P_Vert', 'P_Base' },'RowNames',{'U Bound','Nu Bound'});
optimal_tab

load('LDC_design/Nom_errors_all')
% load('LDC_design/nominal_all_qoi_2')
% load('LDC_design/nominal_all_qoi_2_n8')

1; 

nom_bound = error_bound_vec; nom_bi = err_Bi_vec; nom_low = err_low_vec; 

%%% Print out the bound and bi for the nominal and optimal

% Should do the actual bi-fidelity found at the optimal point. 
% min_bi - result if its the minimum of the bi-fidelity data. 
% min_bi_opt - bi found from min locatin of bound. 
% results_mat = [nom_low, nom_bound, nom_bi, min_bound, min_bi]'*100; 
opt_low = [0, 0, 0, 0, 0]'; 

results_mat = [nom_low*100, nom_bound*100, nom_bi*100, opt_low, min_bound, min_bi_opt, min_bi]'; 
1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot QoI realization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1; 

load('u_meshes/u_64_f_2.mat')

% I could plot mean and variance? 
% or most improved. 

% 200 samples. 65 points
Uf_u = u_matrix_0'; 

% best u mid is at u 1.5579, -0.11579
% best p base is at 0.7157, 2.805

Uc_nom_u = load('LDC_design/Nom_u_mid.mat', 'Uc','Ub','sb');
Ub_nom_u = Uc_nom_u.Ub; 
sb_nom_u = Uc_nom_u.sb; 
Uc_nom_u = Uc_nom_u.Uc; 

Uc_opt_u = load('LDC_design/Opt_u_mid.mat', 'Uc','Ub','sb');
Ub_opt_u = Uc_opt_u.Ub; 
sb_opt_u = Uc_opt_u.sb; 
Uc_opt_u = Uc_opt_u.Uc; 

e_low_u = norm(Uc_opt_u - Uf_u)/norm(Uf_u); 

error_b_nom_u = vecnorm(Ub_nom_u-Uf_u)./vecnorm(Uf_u);
error_b_opt_u = vecnorm(Ub_opt_u-Uf_u)./vecnorm(Uf_u);

[~, index_max_u] = max(error_b_nom_u - error_b_opt_u); 
% norm(Ub_opt_u(:,index_max_u) - Uf_u(:,index_max_u))
% norm(Ub_nom_u(:,index_max_u) - Uf_u(:,index_max_u))

load 'x_64.mat'

x_highfidelity = x_64(:,1); 

plot_index = index_max_u; 

% % plot H, Nominal LF, Nominal BF

figure 
p1 = plot(x_highfidelity,Uf_u(:,plot_index),'color',c1,'LineWidth',LW);
hold on
p2 = plot(x_highfidelity,Uc_nom_u(:,plot_index),':','color',c2,'LineWidth',LW);
p3 = plot(x_highfidelity,Ub_nom_u(:,plot_index),'-.','color',c3,'LineWidth',LW);
p4 = plot(x_highfidelity,Uc_opt_u(:,plot_index),':','color',c4,'LineWidth',LW);
p5 = plot(x_highfidelity,Ub_opt_u(:,plot_index),'-.','color',c5,'LineWidth',LW);
hold off
xlabel('$x$','interpreter','latex','Fontsize',FS)
ylabel('$U_{\mathrm{mid}}$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend([p1,p2,p3,p4,p5],{'HF','Nominal LF','Nominal BF', 'Optimal LF', 'Optimal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3,p4,p5],{'HF','Nominal LF','Nominal BF', 'Optimal LF', 'Optimal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
%title('Mid Velocity Realization','Interpreter','latex')
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'plots/LDC_U_mid_realization','epsc')
end

Uf_pb = u_matrix_4'; 

Uc_nom_pb = load('LDC_design/Nom_p_base.mat', 'Uc','Ub','sb');
Ub_nom_pb = Uc_nom_pb.Ub; 
sb_nom_pb = Uc_nom_pb.sb; 
Uc_nom_pb = Uc_nom_pb.Uc; 

Uc_opt_pb = load('LDC_design/Opt_p_base.mat', 'Uc','Ub','sb');
Ub_opt_pb = Uc_opt_pb.Ub; 
sb_opt_pb = Uc_opt_pb.sb; 
Uc_opt_pb = Uc_opt_pb.Uc; 

e_low_pb = norm(Uc_opt_pb - Uf_pb)/norm(Uf_pb); 

error_b_nom_pb = vecnorm(Ub_nom_pb-Uf_pb)./vecnorm(Uf_pb);
error_b_opt_pb = vecnorm(Ub_opt_pb-Uf_pb)./vecnorm(Uf_pb);

[~, index_max_pb] = max(error_b_nom_pb - error_b_opt_pb); 
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
ylabel('$P_{\mathrm{base}}$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend([p1,p2,p3,p4,p5],{'HF','Nominal LF','Nominal BF', 'Optimal LF', 'Optimal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3,p4,p5],{'HF','Nominal LF','Nominal BF', 'Optimal LF', 'Optimal BF'},'interpreter', 'latex', 'fontsize', FS_leg)
%title('Base Pressure Realization','Interpreter','latex')
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'plots/LDC_P_base_realization','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add optimal low for u mid and p base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1; 

results_mat(4,1) = e_low_u*100; results_mat(4,5) = e_low_pb*100;

results_tab = array2table(results_mat,...
    'VariableNames',{'U_Mid', 'U_Vert', 'P_Mid', 'P_Vert', 'P_Base' },...
    'RowNames',{'Nom LF','Nom Bound','Nom BF', 'Opt LF' , 'Opt Bound', 'Opt BF','Best BF'});
results_tab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot histogram of errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The optimal should be found from the bound. 

%%% Vizualize  data
n_hist = 20; 

figure
hold on
l1 = xline(nom_bound(1)*100,'LineWidth', LW, 'Color', c1);
l2 = xline(min_bound(1),'--','LineWidth', LW, 'Color', c2);

h1 = histogram(abs(error_b_nom_u)*100,n_hist,'FaceColor',c1);
h2 = histogram(abs(error_b_opt_u)*100,n_hist,'FaceColor',c2);
hold off
legend([l1, l2, h1,h2],{'Nominal Bound','Optimal Bound','Nominal Ensemble','Optimal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
yl = ylim;
ylim([yl(1),yl(2)*(1+0.05)]);
xl = xlim;
xlim([xl(1),xl(2)+xl(2)*0.05]);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_1)
%title('U Mid','Interpreter','latex')

if save_on ==1
    saveas(gcf,'plots/LDC_U_mid_hist','epsc')
%     saveas(gcf,'plots/LDC_U_mid_hist','png')
end

figure
hold on
l1 = xline(nom_bound(5)*100,'LineWidth', LW, 'Color', c1);
l2 = xline(min_bound(5),'--','LineWidth', LW, 'Color', c2);

h1 = histogram(abs(error_b_nom_pb)*100,n_hist,'FaceColor',c1);
h2 = histogram(abs(error_b_opt_pb)*100,n_hist,'FaceColor',c2);
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
%title('P Base','Interpreter','latex')

if save_on ==1
    saveas(gcf,'plots/LDC_P_base_hist','epsc')
%     saveas(gcf,'plots/LDC_P_base_hist','png')

end


% Do MID: 
% Uc_nom_u
% sb_nom_u
% 
% Uc_nom_pb
% sb_nom_pb

% r = 1 - 4.1 % (can't go lower)
B = Uc_nom_u/norm(Uc_nom_u,'fro');
[P_s,ix] = matrixIDvR(B,r);
err_Bhat = norm(B-B(:,ix)*P_s) 

% r = 2 %is < 1 % 
B = Uc_nom_pb/norm(Uc_nom_pb,'fro');
[P_s,ix] = matrixIDvR(B,r);
err_Bhat = norm(B-B(:,ix)*P_s) 
    
