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
%%%  nu test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu_orig = 0.01; 
delta_nu_0 = 0; 
delta_nu_1 = -0.5; 

nu_0 = nu_orig*(1+delta_nu_0); 
nu_1 = nu_orig*(1+delta_nu_1); 

x_vec = linspace(0,1,100);
[X,Y] = meshgrid(x_vec); 

%%% nu linear 
nu_mat_1 = zeros(length(x_vec), length(x_vec)); 

for i_x = 1:length(x_vec)
    for i_y = 1:length(x_vec)
        nu_mat_1(i_x,i_y) = nu_0+(nu_1-nu_0)*x_vec(i_y);
    end
end

figure
[c,h]=contourf(X,Y,nu_mat_1',10);
set(h, 'edgecolor','none');
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
pbaspect([1 1 1])
colorbar 
title('Linear $\nu$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots_nu/nu_linear_cont','png')
end

figure
plot([0,1], [nu_0,nu_1],'color',c1,'LineWidth',LW)
xlabel('$y$','interpreter','latex','Fontsize',FS)
ylabel('$\Delta \nu$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
title('Linear','Interpreter','latex')
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'plots_nu/nu_linear_line','png')
end

%%% nu sigmoid 
delta_nu_0 = 1; 
delta_nu_1 = -0.5; 

nu_0 = nu_orig*(1+delta_nu_0); 
nu_1 = nu_orig*(1+delta_nu_1); 

nu_mat_2 = zeros(length(x_vec), length(x_vec)); 

x_vec = linspace(0,1,100);
s_tight = 30; 
s_center = 0.25; 

s_1 = nu_0+(nu_1-nu_0)./(1+exp(s_tight*(x_vec-s_center)));

% vort = [0.6215,0.7357]; 
vort =  [0.6684, 0.7357];

for i_x = 1:length(x_vec)
    for i_y = 1:length(x_vec)
        d_1 = sqrt((x_vec(i_x) - vort(1))^2+(x_vec(i_y) - vort(2))^2);
        
        nu_mat_2(i_x,i_y) = nu_0+(nu_1-nu_0)./(1+exp(s_tight*(d_1-s_center)));
    end
end

figure
[c,h]=contourf(X,Y,nu_mat_2',10);
set(h, 'edgecolor','none');
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
pbaspect([1 1 1])
colorbar 
title('Sigmoid $\nu$','Interpreter', 'latex')

if save_on ==1
    saveas(gcf,'plots_nu/nu_sigmoid_cont','png')
end

% I have nu_0 and nu_1 now do a sigmoid from 1 to the other. Magnitude is
% difference. 

x_vec = linspace(0,1,100);
s_tight = 30; 
s_center = 0.25; 

s_1 = nu_0+(nu_1-nu_0)./(1+exp(s_tight*(x_vec-s_center)));

figure
plot(x_vec, s_1,'color',c1,'LineWidth',LW)
xlabel('distance','interpreter','latex','Fontsize',FS)
ylabel('$\nu$','interpreter','latex','Fontsize',FS)
axis tight
% xlim([0,1]); ylim([0,1]);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
title('Sigmoid','Interpreter','latex')
set(gcf,'Position',size_1)

if save_on ==1
    saveas(gcf,'plots_nu/nu_sigmoid_line','png')
end

% %%% nu sigmoid ring 
% 
% % Plan is to do a picewise sigmoid... have a peak 
% 
% 
% % Could do a gaussian instead? 
% x_vec = 0:0.01:1;
% s_f = 0.1; 
% mu_f = 0.2; 
% 
% f = 1/(s_f.*sqrt(2.*pi)).*exp(-0.5.*((x_vec-mu_f)/s_f).^2);
% 
% 
% figure
% plot(x_vec, f,'color',c1,'LineWidth',LW)
% xlabel('distance','interpreter','latex','Fontsize',FS)
% ylabel('$\nu$','interpreter','latex','Fontsize',FS)
% axis tight
% % xlim([0,1]); ylim([0,1]);
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% title('Offset Gaussian','Interpreter','latex')
% set(gcf,'Position',size_1)



% next try a ring... 
