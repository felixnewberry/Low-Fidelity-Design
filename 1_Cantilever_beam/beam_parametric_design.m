% House keeping
close all
clear all
clc

% Objective:

% Minimize epsilon tau over low-fidelity parameterization of model. 

% Then apply to both lifting technique and basis reduction and compare
% performance. 

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
%%% Choose optimization method and mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% type of search 
individual_search = 0; 
grid_search = 0; % 20x20 - 400 samples
nominal_optimal_tip = 1; 

% Choose parameters to vary
% (mode 0-4 for individual_search, 5 for grid)
% mode = 0; % w
% mode = 1; % h1
% mode = 2; % h2
% mode = 3; % h3
% mode = 4; % h1 = h2
mode = 5; % h1=h2 and h3

% choose delta range
delta_test = 0; % used for individual search of h3 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordinates
load('Beam_data/x_highfidelity.txt')
%Uf fine
load('Beam_data/Uf')
% Uc course
load('Beam_data/Uc')
% xi 
load('Beam_data/xi')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% epsilon tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select subset R columns of N samples from A and B. 
% rank of Uc is r=1. 

r = 1; 
n = r+4; 

nsim = 100; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% individual search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if individual_search == 1
    
% nominal and optimal delta values
% delta_t1_rand = [0, 1.5];
% delta_t3_rand = [0,29]; 

delta_vec = -0.95:0.1:0.95;

if mode == 0
    plot_label = '$ \Delta w [\%]$';
    save_label = 'individual_w_pm95'; 
elseif mode == 1
    plot_label = '$ \Delta h_1 [\%]$';
    save_label = 'individual_h1_pm95'; 
elseif mode == 2
    plot_label = '$ \Delta h_2 [\%]$';
    save_label = 'individual_h2_pm95'; 
elseif mode == 3
    plot_label = '$ \Delta h_3 [\%]$';
    save_label = 'individual_h3_pm95'; 
    if delta_test == 0
        save_label = 'individual_h3_p35';
        delta_vec = 0:1:40;       
    end
elseif mode == 4
    plot_label = '$ \Delta h_1 = \Delta h_2 [\%]$';
    save_label = 'individual_h1h2_p2m1'; 
    delta_vec = -1:0.075:2;
end

% delta_vec = 0; 

error_bound_mat = zeros(length(delta_vec),1); 
error_Bi_mat = zeros(length(delta_vec),1);
efficacy_mat = zeros(length(delta_vec),1);

for i_test = 1:length(delta_vec)
    
[error_bound,err_Bi,efficacy] =  my_beam_bound(delta_vec(i_test),nsim, n, r, mode);

error_bound_mat(i_test) = error_bound;
error_Bi_mat(i_test) =  err_Bi;
efficacy_mat(i_test) = efficacy; 
end

save(strcat('Beam_design/', save_label), 'error_bound_mat', 'delta_vec')

% Plot error bound 
figure
hold on
p1 = plot(100*delta_vec,100*error_bound_mat,'ob-', 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
grid on

% Plot bi-fidelity error 
figure
hold on
p1 = plot(100*delta_vec,100*error_Bi_mat,'sr-', 'LineWidth',LW); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bi $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
grid on

end

if grid_search == 1
%%

mode = 5; 

% 20 x 20 takes 36 s
delta_t3_vec = 0:2:40; 
delta_t1_vec = -1:0.15:2; 

% % 40 x 40 takes 
% delta_t3_vec = 0:1:40; 
% delta_t1_vec = -1:0.075:2; 

tic
error_bound_mat = zeros(length(delta_t1_vec),length(delta_t3_vec)); 
error_Bi_mat = zeros(length(delta_t1_vec),length(delta_t3_vec)); 
efficacy_mat = zeros(length(delta_t1_vec),length(delta_t3_vec)); 

for i_test = 1:length(delta_t1_vec)
    
 fprintf('Grid search completetion: %.2f %% \n', 100*i_test/length(delta_t1_vec)) 
 
for i_t3 = 1:length(delta_t3_vec)

[error_bound,err_Bi,efficacy] =  my_beam_bound([delta_t1_vec(i_test),delta_t3_vec(i_t3)],nsim, n, r, mode);

% results(i_t3,:) = results(i_t3,:) +[error_bound, err_Bi, efficacy]; 
error_bound_mat(i_test,i_t3) = error_bound; 
error_Bi_mat(i_test,i_t3) =err_Bi; 
efficacy_mat(i_test,i_t3) = efficacy; 

% error_efficacy_sum
end

end

save('Beam_design/grid_search','delta_t3_vec','delta_t1_vec',...
    'error_bound_mat', 'error_Bi_mat')

figure
hold on
contourf(100*delta_t3_vec,100*delta_t1_vec,100*error_bound_mat)
colorbar
% 5, 0.2
p1 = plot(0,0,'ro','MarkerSize',8,'linewidth',LW);
% p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',8,'linewidth',LW);
% p3 = plot(100*t3_Bi,100*t1_Bi,'rs','MarkerSize',8,'linewidth',LW);
hold off
% legend([p1,p2,p3],{'Nominal','Optimal', 'Bi'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthWest')
xlabel('$\Delta h_3$ [\%]','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1 = \Delta h_2$ [\%]','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
title('Error bound','Interpreter','latex')


figure
hold on
contourf(100*delta_t3_vec,100*delta_t1_vec,100*error_Bi_mat)
colorbar
p1 = plot(0,0,'ro','MarkerSize',8,'linewidth',LW);
% p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',8,'linewidth',LW);
% p3 = plot(100*t3_Bi,100*t1_Bi,'rs','MarkerSize',8,'linewidth',LW);
hold off
% legend([p1,p2,p3],{'Nominal','Optimal', 'Bi'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthWest')
xlabel('$\Delta h_3 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1= \Delta h_2$ [\%]','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
title('Bi-fidelity error','Interpreter','latex')

end

if nominal_optimal_tip == 1

mode = 6; 
% mode =6 indicates mode 5, and save things. 

% Find nominal and optimal data
delta_t1_vec = [0, 2.0]; 
delta_t3_vec = [0, 40]; 

error_bound_mat = zeros(length(delta_t1_vec),1); 
error_Bi_mat = zeros(length(delta_t1_vec),1); 
efficacy_mat = zeros(length(delta_t1_vec),1); 

for i_test = 1:length(delta_t1_vec)
     
[error_bound,err_Bi,efficacy] =  my_beam_bound([delta_t1_vec(i_test),delta_t3_vec(i_test)],nsim, n, r, mode);

error_bound_mat(i_test) = error_bound; 
error_Bi_mat(i_test) =err_Bi; 
efficacy_mat(i_test) = efficacy; 
end

error_bound_mat
error_Bi_mat
end

% nominal t3 = 5, t2 = t1 = 0.2

