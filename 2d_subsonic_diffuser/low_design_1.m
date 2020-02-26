clear all
close all
clc

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

save_on = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process low-fid AIP QoI - to be moved to simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% New QoI
% vL = importdata('2d_diffuser_data/low_data/v_aip.txt');
% uL = importdata('2d_diffuser_data/low_data/u_aip.txt');
% pL = importdata('2d_diffuser_data/low_data/p_aip.txt');
% probe = importdata('2d_diffuser_data/low_data/p_probe.txt');
% prL = pL./probe; 
% 
% % Plot data
% 
% figure
% plot(vL)
% 
% figure
% plot(uL)
% 
% figure
% plot(pL)
% 
% figure
% plot(prL)
% 
% % Save low-fid_data - would be better if this was pre-computed in sample. 
% % save('2d_diffuser_data/low_nom_raw','vL','uL','pL','prL','probe')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process high-fid AIP QoI - to be moved to simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% New QoI
% vH = importdata('2d_diffuser_data/high_data/v_aip.txt');
% uH = importdata('2d_diffuser_data/high_data/u_aip.txt');
% pH = importdata('2d_diffuser_data/high_data/p_aip.txt');
% probe = importdata('2d_diffuser_data/high_data/p_probe.txt');
% prH = pH./probe; 
% 
% % Plot data
% 
% figure
% plot(vH)
% 
% figure
% plot(uH)
% 
% figure
% plot(pH)
% 
% figure
% plot(prH)
% 
% % Save low-fid_data - would be better if this was pre-computed in sample. 
% % save('2d_diffuser_data/high_raw','vH','uH','pH','prH','probe')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process AIP QoI low and high coordinates. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% h_data = importdata('pvdatata_500_aip_h_example.csv');
% % Check colums 6,7 and 8.
% % x is constant at 0.217342, z is 0. 
% y_h = h_data.data(:,7);

% l_data = importdata('pvdatata_500_aip_l_example.csv');
% Check colums 9, 10, 11 
% y_l = l_data.data(:,10);


% save('2d_diffuser_data/coords_orig','y_l','y_h')

% % Standardize from 0 to 1 

% y_l  = normalize(y_l,'range');
% y_h  = normalize(y_h,'range');
% save('2d_diffuser_data/coords_standardized','y_l','y_h')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Re_order low, high and coords
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('2d_diffuser_data/coords_standardized')
% 
% load('2d_diffuser_data/high_raw')
% load('2d_diffuser_data/low_nom_raw')
% 
% % first order y_h.. 
% [yh, index_h] = sort(y_h);
% [yl, index_l] = sort(y_l);
% 
% % Illustrate distrubution of points with histogram
% 
% n_hist = 20; 
% 
% figure
% hold on
% h1 = histogram(yh,n_hist,'FaceColor',c1);
% h2 = histogram(yl,n_hist,'FaceColor',c2);
% hold off
% legend([h1,h2],{'high','low'},'interpreter', 'latex', 'fontsize', FS_leg)
% xlabel('Standardized $y$','interpreter','latex','Fontsize',FS)
% ylabel('Frequency','interpreter','latex','Fontsize',FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% set(gcf,'Position',size_1)
% 
% %%% Reorder data
% uH = uH(index_h,:); 
% vH = vH(index_h,:); 
% pH = pH(index_h,:); 
% prH = prH(index_h,:); 
% 
% uL = uL(index_l,:); 
% vL = vL(index_l,:); 
% pL = pL(index_l,:); 
% prL = prL(index_l,:); 
% 
% %%% Interpolate
% % 200 evenly spaced points. 
% y_int = linspace(0,1,201);
% 
% uH = interp1(yh, uH, y_int);
% vH = interp1(yh, vH, y_int);
% pH = interp1(yh, pH, y_int);
% prH = interp1(yh, prH, y_int);
% 
% uL = interp1(yl, uL, y_int,'Linear','extrap');
% vL = interp1(yl, vL, y_int,'Linear','extrap');
% pL = interp1(yl, pL, y_int,'Linear','extrap');
% prL = interp1(yl, prL, y_int,'Linear','extrap');
% 
% save('2d_diffuser_data/high_int','vH','uH','pH','prH','probe','y_int')
% save('2d_diffuser_data/low_nom_int','vL','uL','pL','prL','probe', 'y_int')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('2d_diffuser_data/high_int')
% load('2d_diffuser_data/low_nom_int')
% 
% % Plot u, v, p, pr
% n_samps = 100; 
% 
% figure
% subplot(2,2,1)
% p1 = plot(y_int, uH(:,1:n_samps),'-','color',c1,'LineWidth',LW);
% hold on
% p2 = plot(y_int,uL,'-.','color',c2,'LineWidth',LW);
% xlabel('Standardized $y$','interpreter', 'latex', 'fontsize', FS)
% ylabel('$u$ $ms^{-1}$','interpreter', 'latex', 'fontsize', FS)
% axis tight
% title('$u$','interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
% % ylim([y_min, y_max]); 
% 
% legend([p1(1),p2(1)],{'H','L'},'interpreter', 'latex', 'fontsize', FS_leg)
% 
% subplot(2,2,2)
% p1 = plot(y_int, vH(:,1:n_samps),'-','color',c1,'LineWidth',LW);
% hold on
% p2 = plot(y_int,vL,'-.','color',c2,'LineWidth',LW);
% xlabel('Standardized $y$','interpreter', 'latex', 'fontsize', FS)
% ylabel('$v$ $ms^{-1}$','interpreter', 'latex', 'fontsize', FS)
% axis tight
% title('$v$','interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
% % ylim([y_min, y_max]); 
% 
% legend([p1(1),p2(1)],{'H','L'},'interpreter', 'latex', 'fontsize', FS_leg)
% 
% subplot(2,2,3)
% p1 = plot(y_int, pH(:,1:n_samps),'-','color',c1,'LineWidth',LW);
% hold on
% p2 = plot(y_int,pL,'-.','color',c2,'LineWidth',LW);
% xlabel('Standardized $y$','interpreter', 'latex', 'fontsize', FS)
% ylabel('$P$ $[Pa]$','interpreter', 'latex', 'fontsize', FS)
% axis tight
% title('Pressure','interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
% % ylim([y_min, y_max]); 
% 
% legend([p1(1),p2(1)],{'H','L'},'interpreter', 'latex', 'fontsize', FS_leg)
% 
% subplot(2,2,4)
% p1 = plot(y_int, prH(:,1:n_samps),'-','color',c1,'LineWidth',LW);
% hold on
% p2 = plot(y_int,prL,'-.','color',c2,'LineWidth',LW);
% xlabel('Standardized $y$','interpreter', 'latex', 'fontsize', FS)
% ylabel('Pressure Recovery','interpreter', 'latex', 'fontsize', FS)
% axis tight
% title('Pressure Recovery','interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
% % ylim([y_min, y_max]); 
% 
% legend([p1(1),p2(1)],{'H','L'},'interpreter', 'latex', 'fontsize', FS_leg)
% set(gcf, 'Position', size_large)
% 
% 
% if save_on ==1
%     saveas(gcf,'plots/qoi_ensemble','epsc')
%     saveas(gcf,'plots/qoi_ensemble','png')
% end
% 1; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Nominal bi-fidelity 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% why do I have nans? In which variables? 

load('2d_diffuser_data/high_int')
load('2d_diffuser_data/low_nom_int')

%%% Adjust pressure to not be huge

% Make go from -1 to 1? 
% Or subtract median? 

% Try 0 to 1: 
pH = normalize(pH,'range'); 
prH = normalize(prH,'range'); 
pL = normalize(pL,'range'); 
prL = normalize(prL,'range'); 

% Plot u, v, p, pr
n_samps = 100; 

qoi_vec = {'u','v','p','pr'}; 
% I should really put these into one matrix... to avoid nasty string stuff.
% 

r_vec = 1:20; 
% r_vec = 5; 

err_bound_vec = zeros(length(r_vec),length(qoi_vec)); 
err_bi_vec = zeros(length(r_vec),length(qoi_vec)); 
err_low_vec = zeros(length(r_vec),length(qoi_vec)); 

for i_r = 1:length(r_vec)
    
    r = r_vec(i_r); 
    n = r+10; 
    
for i_qoi = 1:length(qoi_vec)
    Uc_name = strcat(qoi_vec{i_qoi},'L');
    Uf_name = strcat(qoi_vec{i_qoi},'H','(:,1:',num2str(n_samps),')');

    Uc = eval(Uc_name);
    Uf = eval(Uf_name);
    
    
%     r = 5; 
%     n = r+10; 

    % Subset of vectors for bi-fidelity error estimate
    rand_sample = 1:n; 

    B = Uc/norm(Uc,'fro');
    A = Uf/norm(Uf,'fro');

    B_R = B(:,rand_sample);
    A_R = A(:,rand_sample);

    % Obtain column skeleton of P
    [P_s,ix] = matrixIDvR(B,r);

    % Error bound inputs
    normC = norm(P_s);
    sb = svd(B); 
    err_Bhat = norm(B-B(:,ix)*P_s); 
    N = 100;

    % % % Compute epsilon tau... 
    [~, ahat_error_est,~, ~,~] = ...
        mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);

    err_bound_vec(i_r,i_qoi) = ahat_error_est/norm(A);
    err_bi_vec(i_r,i_qoi) = norm(A-A(:,ix)*P_s)/norm(A);
    err_low_vec(i_r,i_qoi) = norm(Uc - Uf)/norm(Uf); 
    
    1; 
    
end
end

% %%% For constant r
% results_mat = 100*[err_low_vec; err_bound_vec; err_bi_vec]';
% results_tab = array2table(results_mat, 'VariableNames', {'L','Bound', 'Bi'}, 'RowNames', {'u', 'v', 'p', 'pr'});

% Hard to believe the low-fid error on pressure is 0.1 
% This issue is mitigated if I scale pressure 0 - 1

% plot function of r bi-fidelity estimate and bound. 


figure
p1 = plot(r_vec, 100*err_bound_vec(:,1),'s-','color',c1,'LineWidth',LW);
hold on
p2 = plot(r_vec,100*err_bi_vec(:,1),'s-.','color',c1,'LineWidth',LW);
p3 = plot(r_vec, 100*err_bound_vec(:,2),'d-','color',c2,'LineWidth',LW);
p4 = plot(r_vec,100*err_bi_vec(:,2),'d-.','color',c2,'LineWidth',LW);
p5 = plot(r_vec, 100*err_bound_vec(:,3),'o-','color',c3,'LineWidth',LW);
p6 = plot(r_vec,100*err_bi_vec(:,3),'o-.','color',c3,'LineWidth',LW);
p7 = plot(r_vec, 100*err_bound_vec(:,4),'x-','color',c4,'LineWidth',LW);
p8 = plot(r_vec,100*err_bi_vec(:,4),'x-.','color',c4,'LineWidth',LW);
xlabel('$r$','interpreter', 'latex', 'fontsize', FS)
ylabel('Error $[\%]$','interpreter', 'latex', 'fontsize', FS)
axis tight
% title('Pressure Recovery','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
ylim([5, 75]); 

legend([p1,p2,p3, p4, p5, p6, p7, p8],{'u Bound','u Bi', 'v Bound','v Bi', 'P Bound','P Bi', 'Pr Bound','Pr Bi'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gcf, 'Position', size_1)

if save_on ==1
    saveas(gcf,'plots/qoi_r_bi','epsc')
    saveas(gcf,'plots/qoi_r_bi','png')
end

