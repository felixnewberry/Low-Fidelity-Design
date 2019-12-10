clear all
close all
clc

% Objective

% Low-fidelity error does not indicate improvement in bi-fidelity error. 
% Rather the bound suggests, L_hat - L and H - TL are more important. 

% Two avenues to look into:
% 1. Look at L for better and worse bounds/bi-fidelity. How does the
% ensemble change? What characteristic is relevant? Distribution? - Compare to H. 
% 2. Can I learn some transfromation from one colomn of L to one column of
% H, then apply this to the rest of L prior to interpolative decomposition?


% Investigating what model physics lacks suggests different physics is
% desirable to better approximate u mid to P base. This is interesting. 

save_on = 0; 

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
LW_axis = 1; 

% Colors
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; 
c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880]; 
c6 = [0.3010, 0.7450, 0.9330]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot histogram of errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load u mid data

load('u_meshes/u_64_f_2.mat')
% 200 samples. 65 points
Uf_u = u_matrix_0'; 

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

%%% Load P base data 

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

%%% load x 
load 'x_64.mat'
x_highfidelity = x_64(:,1); 

%%% Plot errors

%%% Vizualize  data
n_hist = 20; 

% check data on this one... 
figure
hold on
h1 = histogram(abs(100*error_b_nom_u),n_hist,'FaceColor',c1);
h2 = histogram(abs(100*error_b_opt_u),n_hist,'FaceColor',c2);
hold off
legend([h1,h2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
title('U Mid','Interpreter','latex')

if save_on ==1
    saveas(gcf,'plots/LDC_U_mid_hist','epsc')
end

figure
hold on
h1 = histogram(abs(100*error_b_nom_pb),n_hist,'FaceColor',c1);
h2 = histogram(abs(100*error_b_opt_pb),n_hist,'FaceColor',c2);
hold off
legend([h1,h2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
title('P Base','Interpreter','latex')

if save_on ==1
    saveas(gcf,'plots/LDC_P_base_hist','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Compare L_nom, L_opt and H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uc_nom_u
% Uc_opt_u
% Uf_u

% 65x200 
% u mid 

n_hist = 30; 

figure
hold on
[p_h,x_h]=hist(abs(100*Uf_u(:)./norm(Uf_u)),n_hist,'Normalization','pdf','FaceColor',c3);
p1 = plot(x_h,p_h/sum(p_h),'LineWidth', LW,'color',c1);
[p_ln,x_ln]=hist(abs(100*Uc_nom_u(:)./norm(Uc_nom_u)),n_hist,'Normalization','pdf','FaceColor',c1);
p2 = plot(x_ln,p_ln/sum(p_ln),'LineWidth', LW,'color',c2);
[p_lo,x_lo] = hist(abs(100*Uc_opt_u(:)./norm(Uc_opt_u)),n_hist,'Normalization','pdf','FaceColor',c2);
p3 = plot(x_lo,p_lo/sum(p_lo),'LineWidth', LW,'color',c3);
hold off
xlabel('Value','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
legend([p1,p2,p3],{'H','L Nominal','L Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
title('u mid','Interpreter','latex')

% these plots compare all 65 points across all 200 samples. What If I just
% plot the curves? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Look into C_L 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate bound + approximation for each L. 

% Plot coefficients, these are the crux to more effective lifting. 

% Then see whether I can apply some transform to L before bi-fidelity... 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 1; 
n = r+2; 
rand_sample = 1:n; 

Uf = Uf_u; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Uc and Uf
Uc = Uc_nom_u; 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

[error_bound_u_nom,err_Bi_u_nom, P_s_u_nom] = my_bound_bi(n,r, A, B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Uc and Uf
Uc = Uc_opt_u; 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

[error_bound_u_opt,err_Bi_u_opt, P_s_u_opt] = my_bound_bi(n,r, A, B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
p1 = plot(P_s_u_nom,'x','MarkerSize', MS,'color',c1);
p2 = plot(P_s_u_opt,'o','MarkerSize', MS,'color',c2);
hold off
xlabel('index','interpreter','latex','Fontsize',FS)
ylabel('coefficient','interpreter','latex','Fontsize',FS)
axis tight
legend([p1,p2],{'C Nominal','C Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
title('u mid coefficients','Interpreter','latex')

figure
hold on
p1 = plot(P_s_u_nom- P_s_u_opt,'x','MarkerSize', MS,'color',c1);
hold off
xlabel('index','interpreter','latex','Fontsize',FS)
ylabel('coefficient difference','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
title('u mid coefficients','Interpreter','latex')

figure
hold on
[p_ln,x_ln]=hist(P_s_u_nom,n_hist,'Normalization','pdf','FaceColor',c1);
p2 = plot(x_ln,p_ln/sum(p_ln),'LineWidth', LW,'color',c2);
[p_lo,x_lo] = hist(P_s_u_opt,n_hist,'Normalization','pdf','FaceColor',c2);
p3 = plot(x_lo,p_lo/sum(p_lo),':','LineWidth', LW,'color',c3);
hold off
xlabel('Value','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
legend([p2,p3],{'C Nominal','C Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
title('u mid','Interpreter','latex')

% Nothing yet. How are the coefficients changing? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Transform L 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Uc and Uf
Uc = Uc_nom_u; 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);
% I have A_R available. 
% A_R and B_R both 65x3 
% Phi* B_R  = A_R ... 

Phi = A_R* pinv(B_R); 

% Tansform L 

B = Phi*B; 
[error_bound_u_nom2,err_Bi_u_nom2, P_s_u_nom2] = my_bound_bi(n,r, A, B);

% Employing pseudoinverese. 
% Purpose is effectively to learn a transform to map L to be much more
% similar to H prior to applying interpolative decomposition. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Uc and Uf
Uc = Uc_opt_u; 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

% I have A_R available. 
% A_R and B_R both 65x3 
% Phi* B_R  = A_R ... 

Phi = A_R* pinv(B_R); 

% Tansform L 

B = Phi*B; 

[error_bound_u_opt2,err_Bi_u_opt2, P_s_u_opt2] = my_bound_bi(n,r, A, B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


res_mat_u =  [error_bound_u_nom, error_bound_u_nom2,error_bound_u_opt,...
    error_bound_u_opt2; ...
    err_Bi_u_nom, err_Bi_u_nom2, err_Bi_u_opt, err_Bi_u_opt2]; 

res_tab_u = array2table(res_mat_u,'VariableNames',{'Nom','Nom2', 'Opt', 'Opt2'}, 'RowNames',{'Bound','Bi'})

% It out performs the low-fidelity design... 
% Try this with further QoI? 

% It drops the bi-fidelity error by more than 50%  - interesting. The bound
% is also very tight --- possibly no longer valid? Though I think it is
% since it is applied to the transformed data. 

% Can I proove that this improves the bi-fidelity estimate? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure
r = 1; 
n = r+2; 
rand_sample = 1:n; 

% Set Uc and Uf
Uc = Uc_nom_pb; 
Uf = Uf_pb; 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

[error_bound_P_nom,err_Bi_P_nom, P_s_P_nom] = my_bound_bi(n,r, A, B);

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);
% I have A_R available. 
% A_R and B_R both 65x3 
% Phi* B_R  = A_R ... 

Phi = A_R* pinv(B_R); 

B = Phi*B;
% Tansform L 

B = Phi*B; 
[error_bound_P_nom2,err_Bi_P_nom2, P_s_P_nom2] = my_bound_bi(n,r, A, B);

% Set Uc and Uf
Uc = Uc_opt_pb; 
Uf = Uf_pb; 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

[error_bound_P_opt,err_Bi_P_opt, P_s_P_opt] = my_bound_bi(n,r, A, B);

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);
% I have A_R available. 
% A_R and B_R both 65x3 
% Phi* B_R  = A_R ... 

Phi = A_R* pinv(B_R); 

% Tansform L 

B = Phi*B; 
[error_bound_P_opt2,err_Bi_P_opt2, P_s_P_opt2] = my_bound_bi(n,r, A, B);

res_mat_P =  [error_bound_P_nom, error_bound_P_nom2,error_bound_P_opt,...
    error_bound_P_opt2; ...
    err_Bi_P_nom, err_Bi_P_nom2, err_Bi_P_opt, err_Bi_P_opt2]; 

res_tab_P = array2table(res_mat_P,'VariableNames',{'Nom','Nom2', 'Opt', 'Opt2'}, 'RowNames',{'Bound','Bi'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%