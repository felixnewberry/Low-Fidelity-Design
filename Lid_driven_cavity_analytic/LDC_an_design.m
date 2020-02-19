%%% 
% Low-fidelity model design 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Classical high-fid
% Analytical low-fid

load('L_nom.mat')

nx = 64; 
r = 5; 
n = r+10; 

% [error_bound,err_bi,err_low] = my_ldc_an_bound(nx, n, r);

% Step through r 1 to 30

r_vec = 1:5; 
bound_vec = zeros(2,length(r_vec)); 
bi_vec = zeros(2,length(r_vec)); 

delta_Re = 0;  

for i_r = 1:length(r_vec)
    r = r_vec(i_r); 
    n = r+10; 

    [bound_vec(:,i_r),bi_vec(:,i_r),err_low] = my_ldc_an_bound(nx, n, r, delta_Re);
end

y_max = max(bound_vec(:));
y_min = min(bi_vec(:));

figure
subplot(1,2,1)
p1 = plot(r_vec, bound_vec(1,:),'-','color',c1,'LineWidth',LW);
hold on
p2 = plot(r_vec,bi_vec(1,:),'--','color',c2,'LineWidth',LW);
xlabel('r','interpreter', 'latex', 'fontsize', FS)
ylabel('Error','interpreter', 'latex', 'fontsize', FS)
axis tight
title('u mid','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
ylim([y_min, y_max]); 

subplot(1,2,2)
p1 = plot(r_vec, bound_vec(2,:),'-','color',c1,'LineWidth',LW);
hold on
p2 = plot(r_vec,bi_vec(2,:),'--','color',c2,'LineWidth',LW);
xlabel('r','interpreter', 'latex', 'fontsize', FS)
ylabel('Error','interpreter', 'latex', 'fontsize', FS)
axis tight
title('v vert','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
ylim([y_min, y_max]); 

legend([p1,p2],{'Bound','Bi'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gcf, 'Position', size_large)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Explore Re
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx = 64; 
r = 2; 
n = r+10; 

% [error_bound,err_bi,err_low] = my_ldc_an_bound(nx, n, r);

% Step through r 1 to 30

delta_Re_vec = [-0.9:0.05:1];  

bound_vec = zeros(2,length(delta_Re_vec)); 
bi_vec = zeros(2,length(delta_Re_vec)); 

for i_re = 1:length(delta_Re_vec)

    [bound_vec(:,i_re),bi_vec(:,i_re),err_low] = my_ldc_an_bound(nx, n, r, delta_Re_vec(i_re));
end

y_max = max(bound_vec(:));
y_min = min(bi_vec(:));

figure
subplot(1,2,1)
p1 = plot(delta_Re_vec, bound_vec(1,:),'-','color',c1,'LineWidth',LW);
hold on
p2 = plot(delta_Re_vec,bi_vec(1,:),'--','color',c2,'LineWidth',LW);
xlabel('$\Delta$ Re','interpreter', 'latex', 'fontsize', FS)
ylabel('Error','interpreter', 'latex', 'fontsize', FS)
axis tight
title('u mid','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
ylim([y_min, y_max]); 

subplot(1,2,2)
p1 = plot(delta_Re_vec, bound_vec(2,:),'-','color',c1,'LineWidth',LW);
hold on
p2 = plot(delta_Re_vec,bi_vec(2,:),'--','color',c2,'LineWidth',LW);
xlabel('$\Delta$ Re','interpreter', 'latex', 'fontsize', FS)
ylabel('Error','interpreter', 'latex', 'fontsize', FS)
axis tight
title('v vert','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
ylim([y_min, y_max]); 

legend([p1,p2],{'Bound','Bi'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gcf, 'Position', size_large)


% Conclusion: create this model in FEniCS 
