clear all
close all
% clc

% Objective:

% Minimize epsilon tau over low-fidelity parameterization of model. 

% Then apply to both lifting technique and basis reduction and compare
% performance. 

% save the correct data and put everything in percent. 

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

% run QoI 1,2,3,4 - then do the nom opt... - and check values first...

% timed run: low fid 200 samples is 2.83 s (+0.9 if time from matlab rather
% than python) high fid is 115 s. 40 times faster.
% High 200 samples 
% how much of this time is spent interpolating? 
% choose QoI
% QoI = 5; 

% record Bound:bi:low
% 0 is u mid v 
% 1 is u vert u 
% 2 is P mid
% 3 is P base           
% 4 is P vert




% % make data for presentation
% plot_ghia = 0; 

% use the test to tune r
% line to set bounds of search
% random search to construct response surface with pce. 
point_test = 1; 
line_search = 0; % 1 for nu, 2 for u
grid_search = 0; % 
random_search = 0; % ie use PCE

nom_opt = 0; % Save data for nominal and optimal runs - have to edit bound too

% to change between new rvs (20% and 10 as opposed to 10% and 5) - edit
% my_ldc_bound, edit the saved variables in paramter (this) script and edit
% the rvs loaded in run_LDC_ensemble.py

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
LW_axis = 2; 

% Colors
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; 
c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880]; 
c6 = [0.3010, 0.7450, 0.9330]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fprintf('**************** Simulation Settings ****************\n')
tic

% This runs, but asks for sudo command. 
pyFi_sudo = ' sudo python3.6 LDC_nonlinear.py';
%pyFi = ' python3.6 LDC_nonlinear.py';

pyFi_ensemble = 'sudo python3.6 LDC_ensemble.py';

% number of repetitions to generate PC response curve
N_sim_PCE = 1; %100; 
% check whether this should be eliminated. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate samples 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 200 samples

% % n_samples = 200; 
% % % rho = 997; 
% % % mu_av = rho/Re; 

Re = 100; 
nu_av = 1/Re; 
u_top_av = 1; 

% % xi_mat = (rand(n_samples,2)-0.5)*2; 
% % nu_vec = nu_av*(1+0.05*xi(:,1)); 
% % u_top_vec = u_top_av*(1+0.1*xi(:,2)); 

% % save('nu_vec','nu_vec')
% % save('u_top_vec','u_top_vec')
% % save('xi_mat','xi_mat')

% % nx = 32 took 133.21 s to run. 
% % nx = 8 took 7.132 s to run. 
% % nx = 6 took 4.94 s to run. 
% % nx = 4 took 3.88 s to run. 

% % % load('LDC_data/nu_vec.mat')
% % % load('LDC_data/u_top_vec.mat')
% % % load('LDC_data/xi_mat.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data, vanilla approach first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordinates x_32 and y_32
load('LDC_data/x_32.mat');
load('LDC_data/y_32.mat');
load('LDC_data/x_all.mat');


% random inputs, xi_mat, nu_vec, u_top_vec
% load('LDC_data/xi_mat.mat');
load('LDC_data/u_nu_vec.mat');

% % Results
% load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/u_matrix_32.mat');
% Uf = u_matrix'; 
% load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/u_matrix_8.mat');
% u_mat_8 = u_matrix'; 
% load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/u_matrix_6.mat');
% u_mat_6 = u_matrix'; 
% load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/u_matrix_4.mat');
% u_mat_4 = u_matrix'; 

% % check pressure whole field: 
% load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/p_field_matrix_int_4.mat');
% Uc = p_field_matrix'; 
% 
% load('home/felixnewberry/Documents/Research/10_low_fidelity_design/u_meshes/p_field_matrix_int_32.mat');
% Uf = p_field_matrix'; 

% epsilon tau setup and test. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Point Test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if point_test == 1
    nx = 4; 
    r = 1; 
    n = r+2; 
    
%     Uc = u_mat_4 ;
%     Uc = p_field_matrix'; 
    
% breaks when delta u is 0.7158 and delta_nu is -0.7... are my
% plots wrong??! As u increases - breaks further. 

% This is high velocity, low viscosity. ie increase in Re. - Newton solver did not converge. 
% If I use umfpack instead of mumps: same result. 
% maybe this is too extreme... check plots
    % grid is square - 20 x 20. I think I have axes swapped. Check. 

% u then nu
% corner case
%     [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r, 2.4, -0.7, -0.7);
    
% I can run 2.4, 0. But not 2.5. If I increase from 1000 to 2000 newton
% iterations - no change. 
% If I change absolute and relative tolerance from -8 and -7 to -6 and -6
% it still does not converge. 

% example where the bound fails with standard RVs. What is going wrong? 
%     [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r, 1.7263, 0.2737, 0.2737);
%     [err_low*100, error_bound*100, err_bi*100]

%     [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r, 0, 0, 0);

    [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r, 0, 0, -0.5);


%     [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r, 0, -0.1, 0.1);

    % u mid best
%     [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r, 1.5579, 0.2737, 0.2737);
    % p base best
%     [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r, 0.2175, 1.8316, 1.8316);
%     [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r, 0, 0, 0);

    [err_low*100, error_bound*100, err_bi*100]

        
% It would be excellent if It could run at a higher Re. 
% For instance the corner of region investigatged is: u + 2.4 nu - 0.7 
% standard implementation is u +- 10 % and nu +- 5 %. 
% Try make work for +- 20 and +- 10. 
% u max = 1.2*(1+2.4), nu_min = 0.009*(1-0.7)
% Nominal Re = 1/0.01 = 100. 
% Corner case Re = u_max/nu_min = 4.08/0.0027 = 1500 (from 100) pretty
% drastic change. 
% I can run 1.2(1+2)/0.01 = 360.. 
    
    
%     save('LDC_design/nominal_all_qoi', 'error_bound', 'err_bi', 'err_low'); 
%     save('LDC_design/nominal_all_qoi_2', 'error_bound', 'err_bi', 'err_low'); 
%     save('LDC_design/nominal_all_qoi_2_n8', 'error_bound', 'err_bi', 'err_low'); 

    % with r = 1 and n = 3: 
% Nominal settings  L, Bo, Bi
% 0 is u mid v      43.3, 5.31, 4.43
% 1 is u vert u     13.3, 1.64, 1.32
% 2 is P mid        57.5, 11.5, 9.4
% 3 is P base       42.5, 5.5, 4.01    
% 4 is P vert       82.0, 11.8, 6.23

% If I increase the rv magnitude..
% 0 is u mid v      43.6, 10.15, 6.59
% 1 is u vert u     13.8, 3.24, 2.0
% 2 is P mid        57.4, 31.68, 21.6
% 3 is P base       43.5, 13.5, 7.01    
% 4 is P vert       82.8, 27.38, 11.9

% If r = 2, n = 4: - bi-fidelity is <5 % for all QoI. 
% Still applicable if limited number or runs. 

% Explore these two scenarious, and if promising then push RVs even
% further.. 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Line Search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if line_search >= 1
    

    % Available hyper-parameters are varying nu and u_top. 
    % start with variation +- 10 
    
    nx = 4; 
%     r = 1; 
    r = 1; 
    n = r+2; 
        
    
%     delta_u_vec = -0.6:0.02:0.0;
%     delta_nu_vec = 0.0:0.1:3;
    

    
    % broke: %     p1 = semilogy(1:length(s),s,'-ob','LineWidth',LW);
    
    if line_search == 1
%         delta_nu_vec = -0.5:0.5:3;
        delta_nu_vec_0 = -0.7:0.1:3;
        delta_nu_vec_1 = delta_nu_vec_0; 
%         delta_nu_vec = -0.75:0.05:-0.7; % -0.8 - newton solver does not converge. same at -0.75


        delta_u_vec = zeros(length(delta_nu_vec_0), 1); 
    elseif line_search == 2
%         delta_u_vec = -0.6:0.3:0.0;
%         delta_u_vec = -0.8:0.1:1.6;
       delta_u_vec = -0.8:0.1:2.4;
        
        delta_nu_vec_0 = zeros(length(delta_u_vec), 1); 
        delta_nu_vec_1 = delta_nu_vec_0; 
    elseif line_search == 3
%         delta_nu_vec_0 = -0.5:0.1:0.5;
        delta_nu_vec_0 = -0.9:0.1:3;

        delta_nu_vec_1 = zeros(length(delta_nu_vec_0), 1); 
        delta_u_vec = delta_nu_vec_1; 
    elseif line_search == 4
%         delta_nu_vec_0 = -0.5:0.1:0.5;
        delta_nu_vec_1 = -0.9:0.1:3;

        delta_nu_vec_0 = zeros(length(delta_nu_vec_1), 1); 
        delta_u_vec = delta_nu_vec_1; 
    end
    
    n_samps = length(delta_nu_vec_0); 
    
%     save('LDC_design/delta_nu_vec','delta_nu_vec')
%     save('LDC_design/delta_u_vec','delta_u_vec')
%     save('LDC_design/delta_P_mid_nu_vec','delta_nu_vec')
%     save('LDC_design/delta_P_mid_u_vec','delta_u_vec')
%     
    error_bound_mat = zeros(5,n_samps); 
    error_Bi_mat = zeros(5,n_samps);    
    
    
    for i_param = 1:n_samps

        [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r,delta_u_vec(i_param),delta_nu_vec_0(i_param), delta_nu_vec_1(i_param));

        error_bound_mat(:,i_param) = error_bound; 
        error_Bi_mat(:,i_param) = err_bi; 
        
        fprintf("Percent complete: %d \n",i_param/n_samps*100);       
    end
%     
 1; 
 




if line_search == 1
    plot_label = '$ \Delta \nu [\%]$';
    delta_vec = delta_nu_vec_0; 
%     save('LDC_design/line_qoi_nu_1','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_nu_2','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_nu_2_n8','error_bound_mat', 'delta_vec')

elseif line_search == 2
    plot_label = '$ \Delta U [\%]$';
    delta_vec = delta_u_vec; 
%     save('LDC_design/line_qoi_u_1','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_u_2','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_u_2_n8','error_bound_mat', 'delta_vec')
elseif line_search == 3
    plot_label = '$ \Delta \nu_0 [\%]$';
    delta_vec = delta_nu_vec_0; 
%     save('LDC_design/line_qoi_u_1','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_u_2','error_bound_mat', 'delta_vec')
    save('LDC_design/line_qoi_nu_y_n4_2','error_bound_mat', 'delta_vec')
elseif line_search == 4
    plot_label = '$ \Delta \nu_1 [\%]$';
    delta_vec = delta_nu_vec_1; 
%     save('LDC_design/line_qoi_u_1','error_bound_mat', 'delta_vec')
%     save('LDC_design/line_qoi_u_2','error_bound_mat', 'delta_vec')
    save('LDC_design/line_qoi_nu_y1_n4_2','error_bound_mat', 'delta_vec')
end

1; 
% Plot error bound 
figure
hold on
p1 = plot(100*delta_vec,100*error_bound_mat(1,:), 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(100*delta_vec,100*error_bound_mat(2,:),'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_bound_mat(3,:), 'Color',c3,'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_vec,100*error_bound_mat(4,:), 'Color',c4,'LineWidth',LW,'MarkerSize',MS); 
p5 = plot(100*delta_vec,100*error_bound_mat(5,:), 'Color',c5,'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Base', '$P$ Vert'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
grid on


% Plot bi-fidelity error 
figure
hold on
p1 = plot(100*delta_vec,100*error_Bi_mat(1,:), 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = plot(100*delta_vec,100*error_Bi_mat(2,:),'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_Bi_mat(3,:), 'Color',c3,'LineWidth',LW,'MarkerSize',MS); 
p4 = plot(100*delta_vec,100*error_Bi_mat(4,:), 'Color',c4,'LineWidth',LW,'MarkerSize',MS); 
p5 = plot(100*delta_vec,100*error_Bi_mat(5,:), 'Color',c5,'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bi $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
legend([p1,p2,p3,p4,p5],{'$U$ Mid','$U$ Vert','$P$ Mid','$P$ Base', '$P$ Vert'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
grid on


%     s = svd(Uc); s = s/s(1); 
%     
%   % Plot svd: 
%     figure
%     p1 = semilogy(1:length(s),s,'-ob','LineWidth',LW);
%     xlabel('index $i$', 'interpreter', 'latex', 'fontsize', FS)
%     ylabel('Normalized Singular Value', 'interpreter', 'latex', 'fontsize', FS)
%     grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
%     
    % 200 samples r = 10. 
    % err_Ahat nx = 8 = 1.55% % error_low is 2.72 % 1.20821 s
    % err_Ahat nx = 6 = 1.58% % error_low is 5.94 
    % err_Ahat nx = 4 = 1.67% % error_low is 18.58 
    
    % LINE PLOT TAKE AWAYS... 
end

if grid_search == 1
    
    nx = 4; 
    r = 1; 
    n = r+2; 
    
    u_lim = [-0.8,2.4]; 
    nu_lim = [-0.7,3.0];
    
    n_grid = 20; 
    
    delta_nu_vec = linspace(nu_lim(1),nu_lim(2),n_grid); 
    delta_u_vec = linspace(u_lim(1),u_lim(2),n_grid); 
    
    delta_u_vec = zeros(1, length(delta_nu_vec)); 
    delta_nu_vec_0 = delta_nu_vec; 
    delta_nu_vec_1 = delta_nu_vec; 
    
    error_bound_mat = zeros(5,length(delta_nu_vec),length(delta_u_vec)); 
    error_Bi_mat = zeros(5,length(delta_nu_vec),length(delta_u_vec)); 
    
%     for i_nu = 1:length(delta_nu_vec)
%         for i_u = 1:length(delta_u_vec)
            
    for i_nu = 1:length(delta_nu_vec_0)
        for i_u = 1:length(delta_nu_vec_1)
            
%             delta_u_vec(i_u)
%             delta_nu_vec(i_nu)
   
            
            [error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r,delta_u_vec(i_u),delta_nu_vec_0(i_nu),delta_nu_vec_0(i_u));
            error_bound_mat(:,i_nu, i_u) = error_bound; 
            error_Bi_mat(:,i_nu, i_u) = err_bi; 
            
            1; 
            
        end 
        fprintf("Percent complete: %d \n",i_nu/length(delta_nu_vec)*100);  
    end
    
    % need to think of a good check on data for runs that didn't
    % complete... - just have u_matrix deleted after each call to error
    % bound. 
%         save('LDC_design/grid_search_2_n8','error_bound_mat', 'error_Bi_mat', 'delta_u_vec','delta_nu_vec')
%         save('LDC_design/grid_search_test','error_bound_mat', 'error_Bi_mat', 'delta_u_vec','delta_nu_vec')
%         save('LDC_design/grid_search_nu_linear_test','error_bound_mat', 'error_Bi_mat', 'delta_u_vec','delta_nu_vec_0','delta_nu_vec_1')

end

if random_search == 1

    % check error with 100 samples... 
n_samps = 500; 



% LDC inputs
nx = 4; 


nx = 4; 
r = 1; 
n = r+2; 

u_lim = [-0.8,2.4]; 
nu_lim = [-0.7,3.0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE fit to random points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_total = n_samps;

% xi_samples is constant for Uf and Uc
% new set of random samples: 
xi_rand = rand(n_samps,3);

% load test
% load('xi_rand_test')


u_rand = u_lim(1)+ (u_lim(2)-u_lim(1)).*xi_rand(:,1); 
nu_rand_0 = nu_lim(1)+ (nu_lim(2)-nu_lim(1)).*xi_rand(:,2); 
nu_rand_1 = nu_lim(1)+ (nu_lim(2)-nu_lim(1)).*xi_rand(:,3); 

% n_samps = 2; 

results_mat = zeros(n_samps,3,5); 
% error_bound_mat = zeros(n_samps,1); 
% error_Ahat_mat = zeros(n_samps,1);
% efficacy_mat = zeros(n_samps,1);

tic
% 
% % look at weird behavior. - turns out I was using tolerance with r
% limited at 3 and it stepped from 1 to 3. 
% 
% u_rand = [0, -0.3] ;
% nu_rand = [0, 0.3] ;

% save velocity and pressure field for these realizations? 

% save field. 

if nom_opt == 1
    if QoI == 0
        u_rand = [0, 0];
        nu_rand_0 = [0, 3];
    elseif QoI == 1
        u_rand = [0, -0.6];
        nu_rand_0 = [0, 2.745];     
    elseif QoI == 2
        u_rand = [0, -0.6];
        nu_rand_0 = [0, 1.932];  
    elseif QoI == 3 
        u_rand = [0, -0.02];
        nu_rand_0 = [0, 0.03]; 
    elseif QoI == 4
        u_rand = [0, -0.6];
        nu_rand_0 = [0, 3];  
    end
end

% change for nom opt... in robust way 
if nom_opt == 1
%     n_start = 1;
%     n_end = 1; 
    
    n_start = 2;
    n_end = 2; 
else
    n_start = 1;
    n_end = n_samps;
end


for i_t = n_start:n_end
    
% [error_bound, err_bi, efficacy] =  my_ldc_bound(QoI, nx,n, r,u_rand(i_t),nu_rand_0(i_t),nu_rand_0(i_t)); 
[error_bound,err_bi,err_low] = my_ldc_bound(nx,n, r,u_rand(i_t),nu_rand_0(i_t), nu_rand_1(i_t));
% error_bound
% err_Ahat
% 1; 

% error_bound = ahat_error_est/norm(Uf);
% err_Ahat = norm(Uf-Uf(:,ix)*P_s)/norm(Uf);
% efficacy = error_bound/err_Ahat;

results_mat(i_t,:,:) = [error_bound,err_bi,err_low]'; 


fprintf("Percent complete: %d \n",i_t/n_samps*100);   
end

% error_bound_mat
% error_Ahat_mat
% efficacy_mat
xi_rand = xi_rand*2-1;

if nom_opt ~= 1
%     save(strcat('LDC_design/rand_',save_label),'error_bound_mat','error_Ahat_mat','efficacy_mat', 'u_lim','nu_lim','xi_rand');
save(strcat('LDC_design/rand_','nu_linear_2'),'results_mat','u_lim','nu_lim','xi_rand');
end

end
