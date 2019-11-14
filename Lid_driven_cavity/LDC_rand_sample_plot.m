

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
%%% Random search 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('LDC_design/rand_nu_linear')

load('LDC_design/nominal_all_qoi_2')
% err_bi and error_bound nominal values - plot these. 

% nu_lim, u_lim, xi_rand and results_mat. 

% 0 is u mid v 
% 1 is u vert u 
% 2 is P mid
% 3 is P base           
% 4 is P vert

% results_mat is 500x3x5 error_bound, bi and low. 



%%% Plot samples

figure
plot3(xi_rand(:,1), xi_rand(:,2), xi_rand(:,3), 'o','MarkerSize',MyMarkerSize)

% sanatize input by replacing zeros (runs that didn't complete) with NaNs. 
% exclude 0s - these are nans. 
results_mat(results_mat==0) = nan;

n_samp = length(results_mat); 

% order in terms of bound to see predictive measure. Different for each
% qoi. 


%%% Plot the nominal, bound and bi - possibly ordered. 

% Step through each qoi 
qoi_vec = 1:5; 

plot_label = ["$U$ Mid","$U$ Vert", "$P$ Mid", "$P$ Base", "$P$ Vert"]; 
plot_save = ["u_Mid","U_Vert", "P_Mid", "P_Base", "P_Vert"]; 

min_bound = zeros(5,1); 
min_bi_opt = zeros(5,1);  
min_bi = zeros(5,1); 

for i_qoi = 1:length(qoi_vec)
[~,idx]  = sort(results_mat(:,1,i_qoi));
results_sort = results_mat(idx,:,i_qoi); 


min_bound(i_qoi) = results_sort(1,1);
min_bi_opt(i_qoi) = results_sort(1,2);
min_bi(i_qoi) = min(results_sort(:,2));

figure
p1 = plot(100*results_sort(:,1),'x','Color',c1,'MarkerSize',MyMarkerSize);
hold on
p2 = plot(100*results_sort(:,2),'o','Color',c2,'MarkerSize',MyMarkerSize);
p3 = plot(100*error_bound(i_qoi)*ones(n_samp,1),'-','Color',c1,'linewidth',LW);
p4 = plot(100*err_bi(i_qoi)*ones(n_samp,1),'-','Color',c2,'linewidth',LW);
hold off
legend([p3,p4,p1,p2],{'Nom Bound','Nom Bi', 'Samp Bound', 'Samp Bi'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('Sample $i$','interpreter','latex','Fontsize',FS)
ylabel('Error [\%]','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_large)
title(strcat(plot_label(i_qoi), ' Samples'),'Interpreter','latex')


if save_on ==1
    saveas(gcf,strcat('Plots_png/LDC_2_rand_samp_', plot_save(i_qoi),'.png'))

end
end

% Plot all 5 qoi then 
% create table

% Just taking the minimum bound for U mid, P Vert (minimal inprovement), U
% Vert and P mid will work well 
% For P Base it appears that some measure of robustness may be necessary to
% locate the best peforming bound that also has the most trust. 

% Ie, what if I take the 200 best peforming and plot them? Tentative idea.
% Entire method really relies on the bound response surface behaing
% similarly to the bi. In this case it doesn't really. 

nom_bound = error_bound; nom_bi = err_bi; nom_low = err_low; 

results_mat_nu_linear = [nom_low, nom_bound, nom_bi, min_bound, min_bi_opt, min_bi]'*100; 

results_tab_nu_linear = array2table(results_mat_nu_linear,...
    'VariableNames',{'U_Mid', 'U_Vert', 'P_Mid', 'P_Base', 'P_Vert' },'RowNames',{'Nom Low','Nom Bound','Nom Bi', 'Opt Bound', 'Opt Bi','Best Bi'});
results_tab_nu_linear
1; 


