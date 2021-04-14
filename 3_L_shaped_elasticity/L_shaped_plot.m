% Felix Newberry
% Date: 01-27-20

% Plot L-shaped elasticity 

clear all
close all
% clc

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot point sample 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Presently this plots an already saved Ub - not anything new
plot_qoi = 2; 

% qoi: 
% 0 - displacement line
% 1 - displacement field
% 2 - stress field

nsim = 200; 
load('L_data/point_test')

load('L_data/dof_coords_f')
dof_coords_f = dof_coords; 
load('L_data/dof_coords_c')
dof_coords_c = dof_coords; 

if plot_qoi == 0
    
    % qoi_0 
    load('L_data/Idx_f')
    load('L_data/Idx_c')

    load('L_data/Uf_line')
    Uf = U(Idx_f,:); 
    load('L_data/Uc_line')
    Uc = U(Idx_c,:); 

    load('L_data/x_f')
    load('L_data/x_c')
    load('L_data/Ub_line')
    
    load('L_data/Uc_line_2')
        
        % interploate surface solution
    Uc_nom_int = zeros(length(x_f),nsim); 

    for i_int = 1:nsim
        Uc_nom_int(:,i_int) = interp1(x_c,Uc(:,i_int),x_f); 
        1; 

    end
elseif plot_qoi == 1
    load('L_data/Uf_field')
    Uf = U; 
    load('L_data/Uc_field')
    Uc = U; 
    load('L_data/Ub_field')
    Ub = U; 
    Uc_nom_int = zeros(length(dof_coords_f),nsim); 

    for i_int = 1:nsim
        F = scatteredInterpolant(dof_coords_c(:,1),dof_coords_c(:,2), Uc(:,i_int));
        Uc_nom_int(:,i_int) = F(dof_coords_f(:,1),dof_coords_f(:,2));
    end
    
elseif plot_qoi == 2
    load('L_data/Uf_stress')
    Uf = Uf; 
    load('L_data/Uc_stress')
    Uc = Uc; 
    load('L_data/Ub_stress')
    Ub = U; 
    Uc_nom_int = zeros(length(dof_coords_f),nsim); 

    for i_int = 1:nsim
        F = scatteredInterpolant(dof_coords_c(:,1),dof_coords_c(:,2), Uc(:,i_int));
        Uc_nom_int(:,i_int) = F(dof_coords_f(:,1),dof_coords_f(:,2));
    end
end

if plot_qoi == 0 
    figure
    p1 = plot(x_f,Uf(:,end),'-x','color',c1, 'LineWidth', LW, 'MarkerSize', MS);
    hold on
    p2 = plot(x_c,Uc(:,end),'-o','color',c2, 'LineWidth', LW, 'MarkerSize', MS);
    % p3 = plot(x_f,Uc_int(:,end),'-d','color',c3, 'LineWidth', LW, 'MarkerSize', MS);
    p4 = plot(x_f,Ub(:,end),'-s','color',c3, 'LineWidth', LW, 'MarkerSize', MS);
    hold off
    xlabel('y', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Horizontal displacement', 'interpreter', 'latex', 'fontsize', FS)
    % legend([p1,p2,p3,p4],{'HF','LF','L_int','BF'},'interpreter', 'latex', 'fontsize', FS_leg)
    legend([p1,p2,p4],{'HF','LF','BF'},'interpreter', 'latex', 'fontsize', FS_leg)
    grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
 axis tight;
    
    % pdf
    [f_f,x_pdf_f] = ksdensity(Uf(1,:)); 
    [f_c,x_pdf_c] = ksdensity(Uc(1,:)); 
    [f_b,x_pdf_b] = ksdensity(Ub(1,:));

    figure % point 0, 1
    p1 = plot(x_pdf_f,f_f , 'color',c1,'LineWidth',LW);
    hold on
    p2 = plot(x_pdf_c,f_c , 'color',c2,'LineWidth',LW);
    p3 = plot(x_pdf_b,f_b , 'color',c3,'LineWidth',LW);
    hold off
    xlabel('Horizontal displacement, $u(0.0,1.0)$', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('plot of $u(0.0,1.0)$', 'interpreter', 'latex', 'fontsize', FS)
    legend([p1,p2,p3],{'HF','LF','BF'},'interpreter', 'latex', 'fontsize', FS_leg)
    grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
axis tight;
    
    
elseif plot_qoi == 1 || plot_qoi == 2
    figure
    p1 = plot3(dof_coords_f(:,1),dof_coords_f(:,2),Uf(:,end),'x','color',c1);
    hold on
    p2 = plot3(dof_coords_c(:,1),dof_coords_c(:,2),Uc(:,end),'o','color',c2);
    p3 = plot3(dof_coords_f(:,1),dof_coords_f(:,2),Uc_nom_int(:,end),'d','color',c3);
    hold off
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
    zlabel('Horizontal displacement', 'interpreter', 'latex', 'fontsize', FS)
    legend([p1,p2,p3],{'HF','LF','L_int'},'interpreter', 'latex', 'fontsize', FS_leg)
    grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
 axis tight;
    
end


error_L = norm(Uf-Uc_nom_int)/norm(Uf);
error_B = norm(Uf-Ub)/norm(Uf);

fprintf("LF error:  %d \n Error Bound:  %d \n BF error:  %d \n",error_L, error_bound, error_Ahat);

% svd
figure
p1 = semilogy(svd(Uf)/max(svd(Uf)),'-x','color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = semilogy(svd(Uc)/max(svd(Uc)),'-o','color',c2,'LineWidth',LW,'MarkerSize',MS); 
% p3 = semilogy(svd(Ub)/max(svd(Ub)),'-x','color',c1,'LineWidth',LW,'MarkerSize',MS); 
xlabel('index $i$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Normalized singular value', 'interpreter', 'latex', 'fontsize', FS)
legend([p1,p2],{'HF','LF'},'interpreter', 'latex', 'fontsize', FS_leg)
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
 axis tight;
xlim([1,30])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot individual sensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot -50 to 50% for each parameter for all three QoI
% Choose parameter
%%% nu 
load('L_design/individual_qoi_0_mp5_p5.mat')

bound_d_line = zeros(6,length(error_bound_mat)); 
bound_d_field = zeros(6,length(error_bound_mat)); 
bound_s_field = zeros(6,length(error_bound_mat)); 

bound_d_line(1,:) = error_bound_mat(:,1);
bound_d_field(1,:) = error_bound_mat(:,2);
bound_s_field(1,:) = error_bound_mat(:,3);

load('L_design/individual_qoi_1_mp5_p5.mat')
bound_d_line(2,:) = error_bound_mat(:,1);
bound_d_field(2,:) = error_bound_mat(:,2);
bound_s_field(2,:) = error_bound_mat(:,3);

load('L_design/individual_qoi_2_mp5_p5.mat')
bound_d_line(3,:) = error_bound_mat(:,1);
bound_d_field(3,:) = error_bound_mat(:,2);
bound_s_field(3,:) = error_bound_mat(:,3);

load('L_design/individual_qoi_3_mp5_p5.mat')
bound_d_line(4,:) = error_bound_mat(:,1);
bound_d_field(4,:) = error_bound_mat(:,2);
bound_s_field(4,:) = error_bound_mat(:,3);

load('L_design/individual_qoi_4_mp5_p5.mat')
bound_d_line(5,:) = error_bound_mat(:,1);
bound_d_field(5,:) = error_bound_mat(:,2);
bound_s_field(5,:) = error_bound_mat(:,3);

load('L_design/individual_qoi_5_mp5_p5.mat')
bound_d_line(6,:) = error_bound_mat(:,1);
bound_d_field(6,:) = error_bound_mat(:,2);
bound_s_field(6,:) = error_bound_mat(:,3);

plot_label = '$ \Delta  [\%],\: \Delta \beta [100$rad$]$';

% Plot for d_line
figure
hold on
p1 = plot(100*delta_vec,100*bound_d_line(1,:),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_vec,100*bound_d_line(2,:),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*bound_d_line(3,:),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_vec,100*bound_d_line(4,:),'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
p5 = plot(100*delta_vec,100*bound_d_line(5,:),'+-', 'Color',c5, 'LineWidth',LW,'MarkerSize',MS); 
p6 = plot(100*delta_vec,100*bound_d_line(6,:),'v--', 'Color',c6, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([p1,p2,p3,p4,p5,p6],{'$\nu \: [\%]$','$E\: [\%]$', '$\ell\: [\%]$', '$\sigma\: [\%]$', '$\beta\: [$rad$]$', '$q\: [\%]$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on

% grid on
set(gcf,'Position',size_1)
%title('QoI: $d$ Individual Sensitivity','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_all_d_line','epsc')
end

% Plot for d_field
figure
hold on
p1 = plot(100*delta_vec,100*bound_d_field(1,:),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_vec,100*bound_d_field(2,:),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*bound_d_field(3,:),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_vec,100*bound_d_field(4,:),'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
p5 = plot(100*delta_vec,100*bound_d_field(5,:),'+-', 'Color',c5, 'LineWidth',LW,'MarkerSize',MS); 
p6 = plot(100*delta_vec,100*bound_d_field(6,:),'v--', 'Color',c6, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([p1,p2,p3,p4,p5,p6],{'$\nu \: [\%]$','$E\: [\%]$', '$\ell\: [\%]$', '$\sigma\: [\%]$', '$\beta\: [$rad$]$', '$q\: [\%]$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on

% grid on
set(gcf,'Position',size_1)
%title('QoI: $d$ field','Interpreter', 'latex')


if save_on ==1
    saveas(gcf,'plots/L_individual_all_d_field','epsc')
end

% Plot for s_field
figure
hold on
p1 = plot(100*delta_vec,100*bound_s_field(1,:),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_vec,100*bound_s_field(2,:),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*bound_s_field(3,:),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_vec,100*bound_s_field(4,:),'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
p5 = plot(100*delta_vec,100*bound_s_field(5,:),'+-', 'Color',c5, 'LineWidth',LW,'MarkerSize',MS); 
p6 = plot(100*delta_vec,100*bound_s_field(6,:),'v--', 'Color',c6, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([p1,p2,p3,p4,p5,p6],{'$\nu \: [\%]$','$E\: [\%]$', '$\ell\: [\%]$', '$\sigma\: [\%]$', '$\beta\: [$rad$]$', '$q\: [\%]$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
% title('QoI: Stress Field','Interpreter', 'latex')

1; 

if save_on ==1
    saveas(gcf,'plots/L_individual_all_s_field','epsc')
end

% for presentation

figure
hold on
p1 = plot(100*delta_vec,100*bound_s_field(1,:),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_vec,100*bound_s_field(2,:),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
%p3 = plot(100*delta_vec,100*bound_s_field(3,:),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
%p4 = plot(100*delta_vec,100*bound_s_field(4,:),'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
%p5 = plot(100*delta_vec,100*bound_s_field(5,:),'+-', 'Color',c5, 'LineWidth',LW,'MarkerSize',MS); 
%p6 = plot(100*delta_vec,100*bound_s_field(6,:),'v--', 'Color',c6, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1,p2],{'$\nu$','$E$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
ylim([20.2247   31.6110])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('QoI: Stress Field','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_all_s_field_nuE','epsc')
end

figure
hold on
p1 = plot(100*delta_vec,100*bound_s_field(1,:),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_vec,100*bound_s_field(2,:),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
%p3 = plot(100*delta_vec,100*bound_s_field(3,:),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
%p4 = plot(100*delta_vec,100*bound_s_field(4,:),'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
%p5 = plot(100*delta_vec,100*bound_s_field(5,:),'+-', 'Color',c5, 'LineWidth',LW,'MarkerSize',MS); 
%p6 = plot(100*delta_vec,100*bound_s_field(6,:),'v--', 'Color',c6, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1],{'$\nu$'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
ylim([20.2247   31.6110])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('QoI: Stress Field','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_all_s_field_nu','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot individual sensitivity search - limits for grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1; 

% nu, \ell, sigma and beta

%%% nu
load('L_design/individual_qoi_nu_0_p7.mat')

figure
hold on
p1 = plot(100*delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_bound_mat(:,3),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$ \Delta  \nu [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$s$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity of $\nu$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_nu','epsc')
end

figure
hold on
% p1 = plot(100*delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_bound_mat(:,3),'o-', 'Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$ \Delta  \nu [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$s$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity of $\nu$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_nu_s_field','epsc')
end


%%% \ell
load('L_design/individual_qoi_cor_0_2p6.mat')

figure
hold on
p1 = plot(100*delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_bound_mat(:,3),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$ \Delta  \ell [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$s$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity of $\ell$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_corr','epsc')
end

figure
hold on
% p1 = plot(100*delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_bound_mat(:,3),'o-', 'Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$ \Delta  \ell [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$s$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity of $\ell$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_corr_s_field','epsc')
end

%%% sig
load('L_design/individual_qoi_sig_mp8_0.mat')

figure
hold on
p1 = plot(100*delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_bound_mat(:,3),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$ \Delta  \sigma [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$s$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity of $\sigma$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_sigma','epsc')
end

figure
hold on
% p1 = plot(100*delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_bound_mat(:,3),'o-', 'Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$ \Delta  \sigma [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$s$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity of $\sigma$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_sigma_s_field','epsc')
end


%%% beta
load('L_design/individual_qoi_theta_m_0_pi_o4.mat')
% load('L_design/individual_qoi_theta_m_pi_o8_pi_o8.mat')

figure
hold on
p1 = plot(delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(delta_vec,100*error_bound_mat(:,3),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$\Delta \beta [$rad$]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$s$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity of $\beta$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_beta','epsc')
end

figure
hold on
% p1 = plot(100*delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(delta_vec,100*error_bound_mat(:,3),'o-', 'Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$\Delta \beta [$rad$]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$s$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
%title('Individual Sensitivity of $\beta$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/L_individual_beta_s_field','epsc')
end

% % plot bound for all QoI... 
% figure
% hold on
% p1 = plot(100*delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(100*delta_vec,100*error_bound_mat(:,3),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
% hold off
% xlabel(plot_label,'interpreter','latex','Fontsize',FS)
% ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% % legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$\sigma$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);%box on
% % grid on
% set(gcf,'Position',size_large)
% 
% if save_on ==1
% %     saveas(gcf,'plots/LDC_1_nu','epsc')
%     saveas(gcf,'Plots/LDC_2_nu','epsc')
% %     saveas(gcf,'plots/LDC_2_nu_8','epsc')
% end
% 
% % plot bi for all QoI
% figure
% hold on
% p1 = plot(100*delta_vec,100*error_bi_mat(:,1),'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_vec,100*error_bi_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(100*delta_vec,100*error_bi_mat(:,3),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
% hold off
% xlabel(plot_label,'interpreter','latex','Fontsize',FS)
% ylabel('Error BF $[\%]$','interpreter','latex','Fontsize',FS)
% % legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Vert', '$P$ Base'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$\sigma$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);%box on
% % grid on
% set(gcf,'Position',size_large)

1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot random sample 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data: 
% load('L_design/nu_corr_sigma_theta_r6')

% xi_rand = xi_rand*2-1;
% 
% n_samps = 200
% 
% %%% PC fit to random points
% N_total = n_samps; 
% N = 200; 
% % nsim_v = N_total - N; 
% % Desired polynomial order of PCE
% d = 4; 
% p = 15; 
% 
% index_pc = nD_polynomial_array(d,p); 
% 
% P = size(index_pc,1);
% 
% 
% 
% 
% clear psi
% for isim=1:N
% %     piset evaluates a multi dimensional pc basis at xi. (legendre 
% %     appropriate for  uniform RV expansion)
%     crow_ref = piset(xi_rand(isim,:),index_pc);
%     psi(isim,:) = crow_ref(1:P);
% end
% 
% i_qoi = 3; 
% 
% % % Solve with least squares
% % c_ref = psi\error_bound_mat(1:N,i_qoi); 
% % c_ref_A = psi\error_Bi_mat(1:N,i_qoi); 
% 
% opts = spgSetParms('iterations',10000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);
% 
% % % solve with spgl1
% weights = get_matrix_weights(psi);
% Psiweights = psi*weights;
% % % sigma is truncation error of PCE (approximated)
% sigma =  cross_val_sigma(psi,error_bound_mat(1:N,i_qoi));
% c_spg = weights*spg_bpdn(Psiweights,error_bound_mat(1:N,i_qoi),sigma*norm(error_bound_mat(1:N,i_qoi)),opts);
% 
% sigma_A =  cross_val_sigma(psi,error_Bi_mat(1:N,i_qoi));
% c_spg_A = weights*spg_bpdn(Psiweights,error_Bi_mat(1:N,i_qoi),sigma_A*norm(error_Bi_mat(1:N,i_qoi)),opts);
% 
% % %%% Start validation
% 
% % clear psi;
% % 
% % for isim=1:nsim_v
% %     crow = piset(xi_rand(isim+N,:),index_pc);
% %     psi(isim,:) = crow(1:P);
% % end
% % 
% % error_val_ls = norm(error_bound_mat(N+1:end,i_qoi)-psi*c_ref)/norm(error_bound_mat(N+1:end,i_qoi));
% % error_val_spg = norm(error_bound_mat(N+1:end,i_qoi)-psi*c_spg)/norm(error_bound_mat(N+1:end,i_qoi));
% % 
% % error_val_ls_A = norm(error_Bi_mat(N+1:end,i_qoi)-psi*c_ref_A)/norm(error_Bi_mat(N+1:end,i_qoi));
% % error_val_spg_A = norm(error_Bi_mat(N+1:end,i_qoi)-psi*c_spg_A)/norm(error_Bi_mat(N+1:end,i_qoi));
% 
% % % PCE stats: 
% fprintf('PCE Statistics: \n');
% % fprintf('LS Bound: %d, LS BF: %d \n',error_val_ls, error_val_ls_A);
% % fprintf('SPG Bound: %d, SPG BF: %d \n',error_val_spg, error_val_spg_A);
% 
% fprintf('SPG Bound: %d, SPG BF: %d \n',sigma, sigma_A);
% 
% % could change this so that spg uses all 200 samples from start. 
% % using least squares and N = 150 vs 50 for validation
% % p = 7
% % 6.9% 
% % p = 8, 7.35 %  
% 
% % Tune PC - then plot. 
% 
% % Can't get below 60% with PC. :/ 

% Stress 
load('L_design/nu_corr_sigma_theta_r6')

%%% Plot

load('L_design/all_nom')

Uc_nom_d_line = Uc_all{1}; 
Ub_nom_d_line = Ub_all{1}; 
sb_nom_d_line = sb_all{1}; 
nom_bi_d_line = error_Bi_all(1); 
nom_bound_d_line = error_bound_all(1); 

Uc_nom_d_field = Uc_all{2}; 
Ub_nom_d_field = Ub_all{2}; 
sb_nom_d_field = sb_all{2}; 
nom_bi_d_field = error_Bi_all(2); 
nom_bound_d_field = error_bound_all(2); 

Uc_nom_s_field = Uc_all{3}; 
Ub_nom_s_field = Ub_all{3}; 
sb_nom_s_field = sb_all{3}; 
nom_bi_s_field = error_Bi_all(3); 
nom_bound_s_field = error_bound_all(3); 

load('L_design/all_opt')
Uc_opt_s_field = Uc_all{3}; 
Ub_opt_s_field = Ub_all{3}; 
% sb_opt_s_field = sb_all{3}; 

% opt_bi_s_field = error_Bi_all(3); 
% opt_bound_s_field = error_bound_all(3); 

[bound_sort_s, bound_index_s] = sort(error_bound_mat(:,3));
bi_sort_s = error_Bi_mat(bound_index_s,3); 

opt_bi_s_field = bi_sort_s(1); 
opt_bound_s_field = bound_sort_s(1); 

delta_mat_s = delta_mat; 

figure
p1 = plot(100*bound_sort_s,'ob','color', c1, 'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(100*bi_sort_s,'xr', 'color', c2,'LineWidth',LW,'MarkerSize',MS); 
p3 = plot([1,200],[nom_bound_s_field*100, nom_bound_s_field*100],'b-','color', c3, 'LineWidth',LW); 
p4 = plot([1,200],[nom_bi_s_field*100, nom_bi_s_field*100],'r--', 'color', c4,'LineWidth',LW); 
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
    saveas(gcf,'plots/L_rand_samples_s','epsc')
end

%%% best performing 
delta_mat_s(:,bound_index_s(1))
bound_sort_s(1)
bi_sort_s(1)

1; 

load('L_data/Uf_stress')
Uf_s = Uf; 
    
% Uc_nom_int = zeros(length(dof_coords_f),nsim); 
% for i_int = 1:nsim
%     F = scatteredInterpolant(dof_coords_c(:,1),dof_coords_c(:,2), Uc_nom_s_field(:,i_int));
%     Uc_nom_int(:,i_int) = F(dof_coords_f(:,1),dof_coords_f(:,2));
% end
% nom_low = norm(Uc_nom_int - Uf_s)/norm(Uf_s); 
% 
% 
% Uc_opt_int = zeros(length(dof_coords_f),nsim); 
% for i_int = 1:nsim
%     F = scatteredInterpolant(dof_coords_c(:,1),dof_coords_c(:,2), Uc_opt_s_field(:,i_int));
%     Uc_opt_int(:,i_int) = F(dof_coords_f(:,1),dof_coords_f(:,2));
% end
% nom_opt = norm(Uc_opt_int - Uf_s)/norm(Uf_s); 

% This has to be checked, seems curious. 
nom_low = norm(Uc_nom_s_field(:,1:nsim) - Uf_s)/norm(Uf_s); 
nom_opt = norm(Uc_opt_s_field(:,1:nsim) - Uf_s)/norm(Uf_s); 
Uc_opt_int = Uc_opt_s_field(:,1:nsim); 


results_mat = [nom_low, nom_bound_s_field, nom_bi_s_field, nom_opt, opt_bound_s_field, opt_bi_s_field]'*100; 

results_tab = array2table(results_mat,...
    'VariableNames',{'S_field' },'RowNames',{'Nom LF','Nom Bound','Nom BF', 'Opt LF', 'Opt Bound', 'Opt BF'});
results_tab

% Check convergence of high fidelity I suppose? Probably necessary to be
% thourough 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot histogram of errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need to checkout the contribution of singularity to the error. 

% The optimal should be found from the bound. 

error_b_nom_s = vecnorm(Ub_nom_s_field-Uf_s)./vecnorm(Uf_s);
error_b_opt_s = vecnorm(Ub_opt_s_field-Uf_s)./vecnorm(Uf_s);
% norm(Ub_opt_s_field-Uf_s)/norm(Uf_s);
%%% Vizualize  data
n_hist = 20; 

% figure
% hold on
% h1 = histogram(abs(100*error_b_nom_s),n_hist,'FaceColor',c1);
% h2 = histogram(abs(100*error_b_opt_s),n_hist,'FaceColor',c2);
% hold off
% legend([h1,h2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
% xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
% ylabel('Frequency','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% set(gcf,'Position',size_1)
% % title('Stress Field','Interpreter','latex')
% 
% 1; 
% 
% if save_on ==1
%     saveas(gcf,'plots/L_s_field_hist','epsc')
% %     saveas(gcf,'plots/LDC_U_mid_hist','png')
% end


% histogram 

results_mat = [nom_low, nom_bound_s_field, nom_bi_s_field, nom_opt, opt_bound_s_field, opt_bi_s_field]'*100; 

x_lim = [0 21.4]; 
y_lim = [0 33]; 

figure
hold on
l1 = xline(nom_bound_s_field*100,'LineWidth', LW, 'Color', c1);
l2 = xline(opt_bound_s_field*100,'--','LineWidth', LW, 'Color', c2);

h1 = histogram(abs(100*error_b_nom_s),n_hist,'FaceColor',c1);
h2 = histogram(abs(100*error_b_opt_s),n_hist,'FaceColor',c2);

hold off
% legend([l1],{'Nominal Bound'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([l1, l2],{'Nominal Bound','Optimal Bound'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([l1, l2, h1],{'Nominal Bound','Optimal Bound','Nominal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
legend([l1, l2, h1,h2],{'Nominal Bound','Optimal Bound','Nominal Ensemble','Optimal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
xlim(x_lim); ylim(y_lim);
yl = ylim;
ylim([yl(1),yl(2)*(1+0.05)]);
xl = xlim;
xlim([xl(1),xl(2)*(1+0.05)]);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf,'Position',size_1)

x_limits = xlim; 
y_limits = ylim; 

if save_on ==1
    saveas(gcf,'plots/L_s_field_hist','epsc')
%     saveas(gcf,'plots/LDC_U_mid_hist','png')
end

% figure
% hold on
% l1 = xline(nom_bound_s_field*100,'LineWidth', LW, 'Color', c1);
% l2 = xline(opt_bound_s_field*100,'--','LineWidth', LW, 'Color', c2);
% 
% h1 = histogram(abs(100*error_b_nom_s),n_hist,'FaceColor',c1);
% h2 = histogram(abs(100*error_b_opt_s),n_hist,'FaceColor',c2);
% 
% hold off
% % legend([l1],{'Nominal Bound'},'interpreter', 'latex', 'fontsize', FS_leg)
% % legend([l1, l2],{'Nominal Bound','Optimal Bound'},'interpreter', 'latex', 'fontsize', FS_leg)
% % legend([l1, l2, h1],{'Nominal Bound','Optimal Bound','Nominal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
% legend([l1, l2, h1,h2],{'Nominal Bound','Optimal Bound','Nominal Ensemble','Optimal Ensemble'},'interpreter', 'latex', 'fontsize', FS_leg)
% xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
% ylabel('Frequency','interpreter','latex','Fontsize',FS)
% axis tight
% xlim(x_limits); 
% ylim(y_limits); 
% % xlim(x_lim); ylim(y_lim);
% % yl = ylim;
% % ylim([yl(1),yl(2)*(1+0.05)]);
% % xl = xlim;
% % xlim([xl(1),xl(2)*(1+0.05)]);
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% set(gcf,'Position',size_1)
% 
% if save_on ==1
%     saveas(gcf,'plots/L_s_field_hist_3','epsc')
% %     saveas(gcf,'plots/LDC_U_mid_hist','png')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot QoI realization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1; 

load('L_data/dof_coords_f')
x_f = dof_coords; 

% Uf_s
% Uc_nom_int
% Uc_opt_int
% Ub_nom_s_field
% Ub_opt_s_field

%%% Relative Errors

Uf = Uf_s; 

% Entire domain
e_L_nom = norm(Uf - Uc_nom_int)/norm(Uf);
e_L_opt = norm(Uf - Uc_opt_int)/norm(Uf);

e_B_nom = norm(Uf - Ub_nom_s_field)/norm(Uf);
e_B_opt = norm(Uf - Ub_opt_s_field)/norm(Uf);

% local to each point (p)
e_L_nom_p = abs((Uf - Uc_nom_int)./(Uf));
e_L_opt_p = abs((Uf - Uc_opt_int)/(Uf));

e_B_nom_p = abs((Uf - Ub_nom_s_field)./(Uf));
e_B_opt_p = abs((Uf - Ub_opt_s_field)./(Uf));

%%% what coordinate sees the biggest improvement? 
e_L_nom_points = vecnorm(Uf - Uc_nom_int,2,2)./vecnorm(Uf,2,2);
e_B_nom_points = vecnorm(Uf - Ub_nom_s_field,2,2)./vecnorm(Uf,2,2);

e_L_opt_points = vecnorm(Uf - Uc_opt_int,2,2)./vecnorm(Uf,2,2);
e_B_opt_points = vecnorm(Uf - Ub_opt_s_field,2,2)./vecnorm(Uf,2,2);

[~,i_points]= max(e_B_nom_points - e_B_opt_points);
 
% from e_B_nom_points(i_points) to e_B_opt_points(i_points)
% what sample sees the biggest improvement?  (16.5 to 4.4)... 

e_L_nom_samples = vecnorm(Uf - Uc_nom_int)./vecnorm(Uf);
e_B_nom_samples = vecnorm(Uf - Ub_nom_s_field)./vecnorm(Uf);

e_L_opt_samples = vecnorm(Uf - Uc_opt_int)./vecnorm(Uf);
e_B_opt_samples = vecnorm(Uf - Ub_opt_s_field)./vecnorm(Uf);

[~,i_sample]= max(e_B_nom_samples - e_B_opt_samples);
% from e_B_nom_samples(i_sample) to e_B_opt_samples(i_sample) (24.6 to
% 4.18)

% need to be careful I don't select a realization that was used to lift - didn't  
e_B_opt_samples(i_sample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pdf 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % pdf of 200 samples at one point... 
% i_points = 1; 
[f_f,x_pdf_f] = ksdensity(Uf(i_points,:)); 
[f_c_nom,x_pdf_c_nom] = ksdensity(Uc_nom_int(i_points,:)); 
[f_b_nom,x_pdf_b_nom] = ksdensity(Ub_nom_s_field(i_points,:));
[f_c_opt,x_pdf_c_opt] = ksdensity(Uc_opt_int(i_points,:)); 
[f_b_opt,x_pdf_b_opt] = ksdensity(Ub_opt_s_field(i_points,:));

figure 
p1 = plot(x_pdf_f,f_f , 'color',c1,'LineWidth',LW);
hold on
p2 = plot(x_pdf_c_nom,f_c_nom , 'color',c2,'LineWidth',LW);
p3 = plot(x_pdf_b_nom,f_b_nom , 'color',c3,'LineWidth',LW);
p4 = plot(x_pdf_c_opt,f_c_opt , 'color',c4,'LineWidth',LW);
p5 = plot(x_pdf_b_opt,f_b_opt , 'color',c5,'LineWidth',LW);
hold off
xlabel('Von Mises Stress at point $(0.0, 0.1903)$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('pdf', 'interpreter', 'latex', 'fontsize', FS)
legend([p1,p2,p3,p4,p5],{'HF','LF Nom','BF Nom','LF Opt','BF Opt'},'interpreter', 'latex', 'fontsize', FS_leg)
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
 axis tight;


% or pdf of all stress for one sample? 
% i_sample = 1; 
[f_f,x_pdf_f] = ksdensity(Uf(:,i_sample)); 
[f_c_nom,x_pdf_c_nom] = ksdensity(Uc_nom_int(:,i_sample)); 
[f_b_nom,x_pdf_b_nom] = ksdensity(Ub_nom_s_field(:,i_sample));
[f_c_opt,x_pdf_c_opt] = ksdensity(Uc_opt_int(:,i_sample)); 
[f_b_opt,x_pdf_b_opt] = ksdensity(Ub_opt_s_field(:,i_sample));

figure 
p1 = plot(x_pdf_f,f_f , 'color',c1,'LineWidth',LW);
hold on
p2 = plot(x_pdf_c_nom,f_c_nom , 'color',c2,'LineWidth',LW);
p3 = plot(x_pdf_b_nom,f_b_nom , 'color',c3,'LineWidth',LW);
p4 = plot(x_pdf_c_opt,f_c_opt , 'color',c4,'LineWidth',LW);
p5 = plot(x_pdf_b_opt,f_b_opt , 'color',c5,'LineWidth',LW);
hold off
xlabel('Von Mises Stress for realization 30', 'interpreter', 'latex', 'fontsize', FS)
ylabel('pdf', 'interpreter', 'latex', 'fontsize', FS)
legend([p1,p2,p3,p4,p5],{'HF','LF Nom','BF Nom','LF Opt','BF Opt'},'interpreter', 'latex', 'fontsize', FS_leg)
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
 axis tight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Possibly I have to do this in paraview? But possibly not. I could inserta
% white box over that corner... 

% can plot actuall stress, or plot the relative errors. 
% Use Uf(:,30), x_f

% why 30?

figure
plot3(x_f(:,1),x_f(:,2),Uf(:,30),'x')

% Interpolate solution - nah. 
[Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);
Uf_surf = griddata(x_f(:,1),x_f(:,2),Uf(:,30),Xq,Yq); 

% to exclude certain points - set to NaN. 
index_1 = find(Xq>0 & Yq>0); 
Uf_surf(index_1) = NaN; 
% V_opt(index_1) = NaN; 

zmax = max(Uf_surf(:)); 


figure
surf(Xq, Yq, Uf_surf,'edgecolor','none')
view(0,90)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
 grid off; box off; axis tight; 
caxis([0,zmax]); 

% Interpolate error nom
[Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);
V_nom = griddata(x_f(:,1),x_f(:,2),e_B_nom_p(:,30),Xq,Yq); 
V_opt = griddata(x_f(:,1),x_f(:,2),e_B_opt_p(:,30),Xq,Yq); 

V_nom(index_1) = NaN; 
V_opt(index_1) = NaN; 

zmax = max(V_nom(:)); 


% figure
% set(gcf, 'Position',  size_2)
% subplot(1,2,1)
% surf(Xq, Yq, V_nom,'edgecolor','none')
% view(0,90)
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
%  grid off; box off; axis tight; 
% caxis([0,zmax]); % could do 0.4
% %title('Nominal Error Realization','Interpreter','latex')
% 
% % zlim([0,zmax])
% 
% % Interpolate error bi
% % [Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);
% 
% subplot(1,2,2)
% surf(Xq, Yq, V_opt,'edgecolor','none')
% view(0,90)
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid off; box off; axis tight;
% caxis([0,zmax]);
% %title('Optimal Error Realization','Interpreter','latex')
% % zlim([0,zmax])
% colorbar
% c =colorbar;
% c.TickLabelInterpreter = 'latex'; 
% set(gcf,'Position',size_2)
% 
% if save_on ==1
%     saveas(gcf,'plots/L_error_Bi_nom_opt','epsc')
% end

figure
set(gcf, 'Position',  size_1)
surf(Xq, Yq, V_nom,'edgecolor','none')
view(0,90)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
 grid off; box off; axis tight; 
 colorbar
c =colorbar;
c.TickLabelInterpreter = 'latex'; 
caxis([0,zmax]); % could do 0.4
%title('Nominal Error Realization','Interpreter','latex')

if save_on ==1
    saveas(gcf,'plots/L_error_Bi_nom','epsc')
end

% zlim([0,zmax])

% Interpolate error bi
% [Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);

figure
set(gcf, 'Position',  size_1)
surf(Xq, Yq, V_opt,'edgecolor','none')
view(0,90)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
grid off; box off; axis tight;
caxis([0,zmax]);
%title('Optimal Error Realization','Interpreter','latex')
% zlim([0,zmax])
colorbar
c =colorbar;
c.TickLabelInterpreter = 'latex'; 

if save_on ==1
    saveas(gcf,'plots/L_error_Bi_opt','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check errors near singularity/contribution to error. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


e_BB_p = e_B_nom_p - e_B_opt_p; % how much worse is nomanl BF compared to Optimal - for all points. 

% 465x200. 
% Plot mean and sd

mu_e_BB_p = mean(e_BB_p,2); 
var_e_BB_p = var(e_BB_p,0,2); 

% Interpolate error nom
[Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);
V_mu = griddata(x_f(:,1),x_f(:,2),mu_e_BB_p,Xq,Yq); 
V_var = griddata(x_f(:,1),x_f(:,2),var_e_BB_p,Xq,Yq); 

% % to exclude certain points - set to NaN. 
V_mu(index_1) = NaN; 
V_var(index_1) = NaN; 

zmax = max(V_mu(:)); 


% figure
% set(gcf, 'Position',  size_2)
% subplot(1,2,1)
% surf(Xq, Yq, V_mu,'edgecolor','none')
% view(0,90)
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid off; box off; axis tight; 
% % caxis([0,zmax]); % could do 0.4
% %title('BF Difference $\mu$','Interpreter','latex')
% 
% % zlim([0,zmax])
% 
% % Interpolate error bi
% % [Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);
% 
% subplot(1,2,2)
% surf(Xq, Yq, V_var,'edgecolor','none')
% view(0,90)
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid off; box off; axis tight;
% % caxis([0,zmax]);
% %title('BF Difference Var','Interpreter','latex')
% % zlim([0,zmax])
% colorbar
% c =colorbar;
% c.TickLabelInterpreter = 'latex'; 
% set(gcf,'Position',size_2)
% 
% if save_on ==1
%     saveas(gcf,'plots/L_error_mu_var','epsc')
% end


figure
set(gcf, 'Position',  size_1)
surf(Xq, Yq, V_mu,'edgecolor','none')
view(0,90)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
grid off; box off; axis tight; 
colorbar
c =colorbar;
c.TickLabelInterpreter = 'latex'; 
caxis([0,zmax]); % could do 0.4


if save_on ==1
    saveas(gcf,'plots/L_error_mu','epsc')
end

% caxis([0,zmax]); % could do 0.4
%title('BF Difference $\mu$','Interpreter','latex')

% zlim([0,zmax])

% Interpolate error bi
% [Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);

figure
set(gcf,'Position',size_1)
surf(Xq, Yq, V_var,'edgecolor','none')
view(0,90)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
grid off; box off; axis tight;
caxis([0,zmax]);
%title('BF Difference Var','Interpreter','latex')
% zlim([0,zmax])
colorbar
c =colorbar;
c.TickLabelInterpreter = 'latex'; 


if save_on ==1
    saveas(gcf,'plots/L_error_var','epsc')
end

Uc_nom_int

B = Uc_nom_int/norm(Uc_nom_int,'fro');
r = 6;
% Gives error above 5% 
[P_s,ix] = matrixIDvR(B,r);
err_Bhat = norm(B-B(:,ix)*P_s) 

