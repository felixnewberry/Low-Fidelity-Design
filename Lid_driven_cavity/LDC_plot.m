% Felix Newberry
% Date: 08-28-19

% Lid driven cavity
% Plot PC results of samples 

clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Task Switches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make data for presentation
plot_ghia = 0; 

plot_mesh = 0; 


%QoI
% choose QoI
QoI = 0; 
% 0 is u_mat
% 1 is P_mat
% 2 is u_vec mid
% 3 is p_vec mid
% 4 is p_vec top

% only works p_top at moment? 
plot_QoI = 0; 
% appropriate for slices ie, mid and top. 

% PCE optimal bound response surface
plot_response = 0; 

plot_p_curious = 0; 
% save bi-fidelity too to compute accurate error and report. Run p_mat now.
% 

save_on = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;     % Line width
MS = 8;     % Marker Size
FS_leg = 16; % Font size legend


size_1 = [0,0,670,515]; 
size_2 = [0,0,1340,515]; 

size_square = [0,0,670,670]; 

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
if plot_ghia == 1
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


% 
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
p2 = plot(u_32, y_32(:,2), 'c-','LineWidth',LW);
p3 = plot(u_8, y_32(:,2), 'b-','LineWidth',LW);
p4 = plot(u_6, y_32(:,2), 'k-','LineWidth',LW);
p5 = plot(u_4, y_32(:,2), 'g-','LineWidth',LW);
p1 = plot(u_ghia_100, y_ghia,'ro-','LineWidth',LW);

xlabel('u velocity', 'interpreter', 'latex', 'fontsize', FS)
xlabel('y', 'interpreter', 'latex', 'fontsize', FS)
%legend('ghia et al', 'fenics')
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
hold off

subplot(2,1,2)
hold on
p2 = plot(x_32(:,1),v_32, 'c-','LineWidth',LW);
p3 = plot(x_32(:,1),v_8, 'b-','LineWidth',LW);
p4 = plot(x_32(:,1),v_6, 'k-','LineWidth',LW);
p5 = plot(x_32(:,1),v_4, 'g-','LineWidth',LW);
p1 = plot(x_ghia, v_ghia_100,'ro-','LineWidth',LW);

xlabel('x', 'interpreter', 'latex', 'fontsize', FS)
ylabel('v velocity', 'interpreter', 'latex', 'fontsize', FS)
legend([p1,p2,p3,p4,p5],{'ghia et al','32','8','6','4'},...
       'interpreter', 'latex', 'fontsize', FS_leg)
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_mesh == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_l = 4; 
n_h = 32; 

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
xlabel('x', 'interpreter', 'latex', 'fontsize', FS)
ylabel('y', 'interpreter', 'latex', 'fontsize', FS)
% grid on; 
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
set(gcf,'Position',size_1)


figure
plot(xx_h, yy_h, 'k','LineWidth',LW/2)
hold on; 
plot(yy_h, xx_h, 'k','LineWidth',LW/2)
xlabel('x', 'interpreter', 'latex', 'fontsize', FS)
ylabel('y', 'interpreter', 'latex', 'fontsize', FS)
% grid on; 
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
set(gcf,'Position',size_1)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot line search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% r = 1, n = 3 ... 
load('LDC_design/line_u_field_nu.mat')
error_u_field_nu = error_bound_mat;
load('LDC_design/line_u_field_u.mat')
error_u_field_u = error_bound_mat;

load('LDC_design/line_P_field_nu.mat')
error_P_field_nu = error_bound_mat;
load('LDC_design/line_P_field_u.mat')
error_P_field_u = error_bound_mat;

load('LDC_design/line_u_mid_nu.mat')
error_u_mid_nu = error_bound_mat;
load('LDC_design/line_u_mid_u.mat')
error_u_mid_u = error_bound_mat;

load('LDC_design/line_p_mid_nu_2.mat')
error_p_mid_nu = error_bound_mat;
load('LDC_design/line_p_mid_u_2.mat')
error_p_mid_u = error_bound_mat;

load('LDC_design/line_p_top_nu.mat')
error_p_top_nu = error_bound_mat;
load('LDC_design/line_p_top_u.mat')
error_p_top_u = error_bound_mat;

load('LDC_design/delta_P_mid_nu_vec.mat')
delta_p_mid_nu_vec = delta_nu_vec; 
load('LDC_design/delta_P_mid_u_vec.mat')
delta_p_mid_u_vec = delta_u_vec; 

load('LDC_design/delta_nu_vec.mat')
load('LDC_design/delta_u_vec.mat')


% I should probably write script so all qoi can be extracted together? 
% maybe... 

% plot nu change for different QoI, ie u_mid, u_field, p_field, p_top
figure
hold on
p1 = plot(100*delta_nu_vec,100*error_u_field_nu,'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_nu_vec,100*error_P_field_nu,'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_nu_vec,100*error_u_mid_nu,'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_nu_vec,100*error_p_top_nu,'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$\Delta \nu [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1,p2,p3,p4],{'$u$ Field','$P$ Field','$u$ Mid','$P$ Top'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
title('Line Search of $\nu$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/LDC_nu_most','epsc')
end


% plot u change for different QoI, ie u_mid, u_field, p_field, p_top
figure
hold on
p1 = plot(100*delta_u_vec,100*error_u_field_u,'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_u_vec,100*error_P_field_u,'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_u_vec,100*error_u_mid_u,'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_u_vec,100*error_p_top_u,'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel('$\Delta U [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1,p2,p3,p4],{'$U$ Field','$P$ Field','$U$ Mid','$P$ Top'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
title('Line Search of $U$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/LDC_u_most','epsc')
end

% plot u change for QoI P Mid
figure
p1 = plot(100*delta_P_mid_u_vec,100*error_p_mid_u,'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
xlabel('$\Delta U [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4],{'$U$ Field','$P$ Field','$U$ Mid','$P$ Top'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
title('Line Search of $U$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/LDC_u_P_mid','epsc')
end

% plot nu change for QoI P Mid
figure
p1 = plot(100*delta_P_mid_nu_vec,100*error_p_mid_nu,'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
xlabel('$\Delta \nu [\%]$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4],{'$U$ Field','$P$ Field','$U$ Mid','$P$ Top'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
title('Line Search of $\nu$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots/LDC_nu_P_mid','epsc')
end

if plot_QoI == 1
    
    if QoI == 4

     pplot_QoI = 1; 
        
    load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/p_line_matrix_int_4.mat');
    Uc = p_line_matrix'; 
    % Uc = [u_field_matrix,p_field_matrix]';
    Uc_orig = Uc; 

    load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/p_line_matrix_top_32.mat');
    Uf = p_line_matrix'; 
    
    load('Uc_trial')
    Uc_opt = Uc_trial; 
            
    else
        fprintf("Unable to plot QoI")
        pplot_QoI = 0;
    end
    
    if pplot_QoI == 1
    % Plot Range to have indication: 
    Uc_mean = mean(Uc_orig,2)';
    Uc_max = max(Uc_orig,[],2)';
    Uc_min = min(Uc_orig,[],2)';
    
    Uf_mean = mean(Uf,2)'; 
    Uf_max = max(Uf,[],2)';
    Uf_min = min(Uf,[],2)';
    Uf_yy = [Uf_min,fliplr(Uf_max)];
    
%     r = 10; nx = 4; 
%     [error_bound,err_Ahat,efficacy] = my_ldc_field_bound(nx,r, -0.9, 5);
    
    
    
    Ucc_mean = mean(Uc_opt,2)'; 
    Ucc_max = max(Uc_opt,[],2)';
    Ucc_min = min(Uc_opt,[],2)';
    Ucc_yy = [Ucc_min,fliplr(Ucc_max)];

    1; 
    load('ensemble_inputs/x_32.mat');
    x_highfidelity = x_32(:,1)'; 
    xx = [x_highfidelity,fliplr(x_highfidelity)];

    figure
    hold on
    p1 = semilogy(x_highfidelity,abs(Uf_mean),'b','LineWidth',LW);
    hold on
    p2 = semilogy(x_highfidelity,abs(Uc_mean),'r','LineWidth',LW);
    p3 = semilogy(x_highfidelity,abs(Ucc_mean),'g','LineWidth',LW);
    set(gca, 'YScale', 'log')
    hold off
    xlabel('$x$','interpreter','latex','Fontsize',FS)
    ylabel('$\vert P_{top} \vert$','interpreter','latex','Fontsize',FS)
    axis tight
    set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
    legend([p1,p2,p3],{'H','Naive L','Optimal L'},'interpreter', 'latex', 'fontsize', FS_leg)
    grid on; 
    
    load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/p_line_matrix_mid_32.mat');
    Uf_mid = p_line_matrix'; 
    Uf_mean_mid = mean(Uf_mid,2)'; 
    
    % use for new QoI stuff. 
%     
%     figure
%     hold on
%     p1 = semilogy(x_highfidelity,abs(Uf_mean_mid),'b','LineWidth',LW);
%     hold on
% %     p2 = semilogy(x_highfidelity,abs(Uc_mean),'r','LineWidth',LW);
% %     p3 = semilogy(x_highfidelity,abs(Ucc_mean),'g','LineWidth',LW);
%     set(gca, 'YScale', 'log')
%     hold off
%     xlabel('$x$','interpreter','latex','Fontsize',FS)
%     ylabel('$\vert P_{mid} \vert$','interpreter','latex','Fontsize',FS)
%     axis tight
%     set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% %     legend([p1,p2,p3],{'H','Naive L','Trial L w'},'interpreter', 'latex', 'fontsize', FS_leg)
%     grid on; 
    end
end

if plot_response == 1 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE fit to random points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if QoI == 0
    load('PC_results_vel_mat')
elseif QoI == 1
    load('PC_results_P_mat')
elseif QoI == 2
%     load('PC_results_vel_midr_10')
    load('PC_results_vel_mid')

elseif QoI == 3
    load('PC_results_P_mid')
elseif QoI == 4
    load('PC_results_P_top')
end

1; 

N_total = length(xi_rand);

N = 150; 
nsim_v = N_total - N; 
% Desired polynomial order of PCE
d = 2; 
p = 7; 

index_pc = nD_polynomial_array(d,p); 

P = size(index_pc,1);

psi = zeros(N,P);

n_reps = 1; 
error_ls = zeros(n_reps,1); 
error_spg = zeros(n_reps,1); 

for i_rep = 1:n_reps

for isim=1:N
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_ref = piset(xi_rand(isim,:),index_pc);
    psi(isim,:) = crow_ref(1:P);
end


% Solve with least squares
c_ref = psi\error_bound_mat(1:N,1); 
c_ref_A = psi\error_Ahat_mat(1:N,1); 

% opts = spgSetParms('iterations',10000,'verbosity',0);
opts = spgSetParms('iterations',10000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);

% % solve with spgl1
weights = get_matrix_weights(psi);
Psiweights = psi*weights;
% % sigma is truncation error of PCE (approximated)
sigma =  cross_val_sigma(psi,error_bound_mat(1:N,1));
c_spg = weights*spg_bpdn(Psiweights,error_bound_mat(1:N,1),sigma*norm(error_bound_mat(1:N,1)),opts);

sigma_A =  cross_val_sigma(psi,error_Ahat_mat(1:N,1));
c_spg_A = weights*spg_bpdn(Psiweights,error_Ahat_mat(1:N,1),sigma*norm(error_Ahat_mat(1:N,1)),opts);


%%% Start validation

clear psi;

for isim=1:nsim_v
    crow = piset(xi_rand(isim+N,:),index_pc);
    psi(isim,:) = crow(1:P);
end

error_val_ls = norm(error_bound_mat(N+1:end,1)-psi*c_ref)/norm(error_bound_mat(N+1:end,1));
error_val_spg = norm(error_bound_mat(N+1:end,1)-psi*c_spg)/norm(error_bound_mat(N+1:end,1));

error_val_ls_A = norm(error_Ahat_mat(N+1:end,1)-psi*c_ref_A)/norm(error_Ahat_mat(N+1:end,1));
error_val_spg_A = norm(error_Ahat_mat(N+1:end,1)-psi*c_spg_A)/norm(error_Ahat_mat(N+1:end,1));


error_ls(i_rep) = error_val_ls
error_spg(i_rep) = error_val_spg

end

% % PCE stats: 
fprintf('PCE Bound Errors: \n');
fprintf('LS Mean: %d, LS Sd: %d \n',mean(error_ls), std(error_ls));
fprintf('SPG Mean: %d SPG Sd: %d \n',mean(error_spg), std(error_spg));  

% make a choive availble for which coefficients to use? for now use LS not
% spg. 

% Surface from PCE

% u_lim = [-0.85,0]; 
% nu_lim = [0,5];

delta_u_vec = u_lim(1):0.01:u_lim(2); 
delta_nu_vec = nu_lim(1):0.05:nu_lim(2);

[X,Y] = meshgrid(delta_u_vec,delta_nu_vec);
Z = zeros(size(X)); 

[xx,yy] = size(X); 
% Transfrom to vector: 
XX = X(:); YY = Y(:); 

% Want -1 to 1 RV. 
%     t1_lim = [0.005,0.6]; 
%     t3_lim = [5,200];

range_u = u_lim(2)-u_lim(1); 
range_nu = nu_lim(2)-nu_lim(1); 

XX_xi = (XX-u_lim(1))/range_u*2-1;
YY_xi = (YY-nu_lim(1))/range_nu*2-1;
xi_grid = [XX_xi,YY_xi]; 

clear psi 

for isim=1:length(XX)
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_ref = piset(xi_grid(isim,:),index_pc);
    psi(isim,:) = crow_ref(1:P);
end

if QoI == 3
    c_use = c_spg; 
    c_use_A = c_spg_A; 
else
    c_use = c_ref; 
    c_use_A = c_ref_A;    
end

% estimate response surface
ZZ = psi*c_use;
ZZ_A = psi*c_use_A;

% Reshape as matrix
Z = reshape(ZZ,[xx,yy]);
Z_A = reshape(ZZ_A,[xx,yy]);



% calculate minimum values and compare to naive
min_bound = min(Z(:));
[i_bound,j_bound] = find(Z==min_bound);

Ahat_from_bound = Z_A(i_bound,j_bound); 

min_Ahat = min(Z_A(:));
[i_Ahat,j_Ahat] = find(Z_A==min_Ahat);

[i_meh,j_naive] = find(X==0);
[i_naive,i_meh] = find(Y==0);
j_naive = j_naive(1); 
i_naive = i_naive(1); 
naive_bound = Z(i_naive,j_naive); 
naive_Ahat = Z_A(i_naive,j_naive); 

% should I report the pc interpolated values of the true values? PC when
% error is low 
pc_res_tab = array2table([naive_bound,naive_Ahat;min_bound,Ahat_from_bound],...
    'VariableNames',{'Bound','Bi'},'RowNames',{'Naive','Optimal'})

fprintf('PC validation error of %d percent for ls \n',error_val_ls*100);
fprintf('PC validation error of %d percent for spg \n',error_val_spg*100);
c_low = min([Z(:);Z_A(:)]); 
c_high = max([Z(:);Z_A(:)]); 

% make colorbar limits consistent (change back if only one plot is used. 
figure
contourf(Y,X,Z)
colorbar
hold on; 
p1 = plot(0,0,'rs','MarkerSize',8,'linewidth',LW);
p2 = plot(delta_nu_vec(i_bound),delta_u_vec(j_bound),'ro','MarkerSize',8,'linewidth',LW);
p3 = plot(delta_nu_vec(i_Ahat),delta_u_vec(j_Ahat),'rx','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Naive','Bound Min', '$\hat{A}$ Min'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('$\Delta \nu$','interpreter','latex','Fontsize',FS)
ylabel('$\Delta u$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
title('Error bound','Interpreter','latex')
% to make comparison to earlier contour
% caxis([c_low,c_high])

% should I make a pce of Ahat too? Yes. 
% validation error for the bound pce was quite excellent, 7.05e-4; 

figure
contourf(Y,X,Z_A)
hold on
colorbar
p1 = plot(0,0,'rs','MarkerSize',8,'linewidth',LW);
p2 = plot(delta_nu_vec(i_bound),delta_u_vec(j_bound),'ro','MarkerSize',8,'linewidth',LW);
p3 = plot(delta_nu_vec(i_Ahat),delta_u_vec(j_Ahat),'rx','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Naive','Bound Min', '$\hat{A}$ Min'},'interpreter', 'latex', 'fontsize', FS_leg)   
xlabel('$\Delta \nu$','interpreter','latex','Fontsize',FS)
ylabel('$\Delta u$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
title('Error A hat','Interpreter','latex')
% caxis([c_low,c_high])




end

if plot_p_curious == 1
    
load('PC_results_P_mid')

load('Uc_p_curious_nom')
Uc_nom = Uc; 

load('Uc_p_curious_better')
Uc_better = Uc; 

Uf= load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/p_line_matrix_mid_32.mat');
Uf = Uf.p_line_matrix';
    Uf= load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/p_line_matrix_mid_32.mat');
Uf = Uf.p_line_matrix';
u_rand = u_lim(1)+ (u_lim(2)-u_lim(1)).*(xi_rand(:,1)+1)/2; 
nu_rand = nu_lim(1)+ (nu_lim(2)-nu_lim(1)).*(xi_rand(:,2)+1)/2; 

%%% trouble shoot error bound directly. 
figure
plot3(nu_rand, u_rand,error_bound_mat,'bo')
hold on
plot3(nu_rand, u_rand,error_Ahat_mat,'rx')
xlabel('$\Delta \nu$','interpreter','latex','Fontsize',FS)
ylabel('$\Delta u$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
title('Error bound','Interpreter','latex')

load('ensemble_inputs/x_32.mat');
x_highfidelity = x_32(:,1)'; 
    
figure
p11 = plot(x_highfidelity,Uf,'b','LineWidth',LW);
hold on
p1 = plot(x_highfidelity,Uf(:,1),'b','LineWidth',LW);

p22 = plot(x_highfidelity,Uc_nom,'r','LineWidth',LW);
p2 = plot(x_highfidelity,Uc_nom(:,1),'r','LineWidth',LW);

p33 = plot(x_highfidelity,Uc_better,'g','LineWidth',LW);
p3 = plot(x_highfidelity,Uc_better(:,1),'g','LineWidth',LW);

% set(gca, 'YScale', 'log')
hold off
xlabel('$x$','interpreter','latex','Fontsize',FS)
ylabel('$\vert P_{mid} \vert$','interpreter','latex','Fontsize',FS)

set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1,p2,p3],{'H','Naive L','Optimal L'},'interpreter', 'latex', 'fontsize', FS_leg)
grid on; 

figure
semilogy(svd(Uc_nom),'rx')
hold on
semilogy(svd(Uc_better),'bo')
set(gca, 'YScale', 'log')
hold on
axis tight

hold off

load('nom_opt/bound_nom.mat') 
true_naive_bound = error_bound; 
load('nom_opt/bound_opt.mat')
true_opt_bound = error_bound; 

load('nom_opt/Ahat_nom.mat') 
true_naive_Ahat = err_Ahat; 
load('nom_opt/Ahat_opt.mat')
true_opt_Ahat = err_Ahat; 

res_tab_data = array2table([true_naive_bound,true_naive_Ahat;true_opt_bound,true_opt_Ahat],...
    'VariableNames',{'Bound','Bi'},'RowNames',{'Naive','Optimal'})
end
    
