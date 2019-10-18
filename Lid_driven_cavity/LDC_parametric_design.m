clear all
close all
clc

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
QoI = 5; 
% 0 is u_mat
% 1 is P_mat
% 2 is u_vec mid
% 3 is p_vec mid
% 4 is p_vec top
% 5 is p_vec mid (vertical)

% % make data for presentation
% plot_ghia = 0; 

% use the test to tune r
% line to set bounds of search
% random search to construct response surface with pce. 
point_test = 0; 
line_search = 0;
random_search = 1; % ie use PCE

nom_opt = 0; % Save data for nominal and optimal runs - have to edit bound too

% QoI 3 needs special treatment, ie change range. 
% Other parameters could explore further along decrease velocity? 

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
load('LDC_data/xi_mat.mat');
load('LDC_data/nu_vec.mat');
load('LDC_data/u_top_vec.mat');

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
    
    [error_bound,err_Ahat,efficacy] = my_ldc_bound(QoI, nx,n, r, 0, 0);
    
%     % 200 samples r = 10. 
%     % err_Ahat nx = 8 = 1.55% % error_low is 2.72 % 1.20821 s
%     % err_Ahat nx = 6 = 1.58% % error_low is 5.94 
%     % err_Ahat nx = 4 = 1.67% % error_low is 18.58 
    
    % bi-fid estimate still
    % good!
    
    error_bound
    err_Ahat
%    	error_low = norm(Uf-Uc)/norm(Uf)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Line Search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if line_search >= 1
    

    % Available hyper-parameters are varying nu and u_top. 
    % start with variation +- 10 
    
    nx = 4; 
    r = 1; 
    n = r+2; 
    
    % r = 1, n = 3 seems okay for U and P fields, U mid and P top. 
    % what's appropriate for P mid? 

    % accidentily ran all samples as n = r+1... 
    % calibrate now. 
    
%     delta_u_top_vec = -0.75:0.02:0.75;
%     delta_nu_vec = -0.1:0.02:0.1;
    % broke: %     p1 = semilogy(1:length(s),s,'-ob','LineWidth',LW);
    
    if line_search == 1
        delta_nu_vec = 0.0:0.1:3;
        
        if QoI == 3 
            delta_nu_vec = [0.0:0.002:0.03];
        end
        
        delta_u_vec = zeros(length(delta_nu_vec), 1); 
    elseif line_search == 2
        delta_u_vec = [-0.6:0.02:0.0];
        
        if QoI == 3 
%             delta_u_vec = [-0.02:0.001:0.0];
            delta_u_vec = [-0.4, -0.02, 0.0];
        end
            
        delta_nu_vec = zeros(length(delta_u_vec), 1); 
    end
    
    
%     save('LDC_design/delta_nu_vec','delta_nu_vec')
%     save('LDC_design/delta_u_vec','delta_u_vec')
%     save('LDC_design/delta_P_mid_nu_vec','delta_nu_vec')
%     save('LDC_design/delta_P_mid_u_vec','delta_u_vec')
    
    error_bound_mat = zeros(1,length(delta_nu_vec)); 
    error_Bi_mat = zeros(1,length(delta_nu_vec));    
    
    
    for i_param = 1:length(delta_nu_vec)

        [error_bound,err_Ahat,efficacy] = my_ldc_bound(QoI, nx,n, r,delta_u_vec(i_param),delta_nu_vec(i_param));

        error_bound_mat(i_param) = error_bound; 
        error_Bi_mat(i_param) = err_Ahat; 
        
        fprintf("Percent complete: %d \n",i_param/length(delta_nu_vec)*100);       
    end
    
1; 




if line_search == 1
    plot_label = '$ \Delta \nu [\%]$';
    delta_vec = delta_nu_vec; 
    if QoI == 0
        save('LDC_design/line_u_field_nu','error_bound_mat')
    elseif QoI == 1
        save('LDC_design/line_P_field_nu','error_bound_mat')
    elseif QoI == 2
        save('LDC_design/line_u_mid_nu','error_bound_mat')
    elseif QoI == 3
        save('LDC_design/line_p_mid_nu_3','error_bound_mat')
    elseif QoI == 4
        save('LDC_design/line_p_top_nu','error_bound_mat')
    elseif QoI == 5
        save('LDC_design/line_p_mid_vert_nu','error_bound_mat')
    end
elseif line_search == 2
    plot_label = '$ \Delta U [\%]$';
    delta_vec = delta_u_vec; 
    if QoI == 0
        save('LDC_design/line_u_field_u','error_bound_mat')
    elseif QoI == 1
        save('LDC_design/line_P_field_u','error_bound_mat')
    elseif QoI == 2
        save('LDC_design/line_u_mid_u','error_bound_mat')
    elseif QoI == 3
        save('LDC_design/line_p_mid_u_3','error_bound_mat')
    elseif QoI == 4
        save('LDC_design/line_p_top_u','error_bound_mat')
    elseif QoI == 5
        save('LDC_design/line_p_mid_vert_u','error_bound_mat')
    end
end

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
end

if random_search == 1

    % check error with 100 samples... 
n_samps = 100; 

% 0 is u_mat
% 1 is P_mat
% 2 is u_vec mid
% 3 is p_vec mid
% 4 is p_vec top

if QoI == 0
    u_lim = [-0.6,0]; 
    nu_lim = [0,3.0];
    r = 1; 
    n = r+2; 
    save_label = 'u_field';
elseif QoI == 1
    u_lim = [-0.6,0]; 
    nu_lim = [0,3.0];
    r = 1; 
    n = r+2; 
    save_label = 'P_field';
elseif QoI == 2
    u_lim = [-0.6,0]; 
    nu_lim = [0,3.0];
    r = 1; 
    n = r+2; 
    save_label = 'u_mid';
elseif QoI == 3
    u_lim = [-0.02,0]; 
    nu_lim = [0,0.03];
    r = 1;   
    n = r+2; 
    save_label = 'P_mid';
elseif QoI == 4
    u_lim = [-0.6,0]; 
    nu_lim = [0,3.0];
    r = 1;    
    n = r+2; 
    save_label = 'P_top';
elseif QoI == 5
    u_lim = [-0.6,0]; 
    nu_lim = [0,3.0];
    r = 1;    
    n = r+2; 
    save_label = 'P_mid_vert';
end

% LDC inputs
nx = 4; 
% r = 10; 

% r = 2; % try this for mid-plane

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE fit to random points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_total = n_samps;

% N = 150; 
% nsim_v = N_total - N; 
% % Desired polynomial order of PCE
% d = 2; 
% p = 7; 

% index_pc = nD_polynomial_array(d,p); 

% P = size(index_pc,1);


% psi = zeros(N,P);

% n_reps = 1; 
% error_ls = zeros(n_reps,1); 
% error_spg = zeros(n_reps,1); 

% Save outputs to file because this is expensive. 



% for i_rep = 1:n_reps

% These samples are to explore u_L and nu_L
% xi_samples is constant for Uf and Uc
% new set of random samples: 
xi_rand = rand(n_samps,2);

% load test
% load('xi_rand_test')


u_rand = u_lim(1)+ (u_lim(2)-u_lim(1)).*xi_rand(:,1); 
nu_rand = nu_lim(1)+ (nu_lim(2)-nu_lim(1)).*xi_rand(:,2); 

% n_samps = 2; 

results = zeros(n_samps,3); 
error_bound_mat = zeros(n_samps,1); 
error_Ahat_mat = zeros(n_samps,1);
efficacy_mat = zeros(n_samps,1);

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
        nu_rand = [0, 3];
    elseif QoI == 1
        u_rand = [0, -0.6];
        nu_rand = [0, 2.745];     
    elseif QoI == 2
        u_rand = [0, -0.6];
        nu_rand = [0, 1.932];  
    elseif QoI == 3 
        u_rand = [0, -0.02];
        nu_rand = [0, 0.03]; 
    elseif QoI == 4
        u_rand = [0, -0.6];
        nu_rand = [0, 3];  
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
    
[error_bound, err_Ahat, efficacy] =  my_ldc_bound(QoI, nx,n, r,u_rand(i_t),nu_rand(i_t)); 

% error_bound
% err_Ahat
% 1; 

% error_bound = ahat_error_est/norm(Uf);
% err_Ahat = norm(Uf-Uf(:,ix)*P_s)/norm(Uf);
% efficacy = error_bound/err_Ahat;

error_bound_mat(i_t) = error_bound; 
error_Ahat_mat(i_t) =  err_Ahat; 
efficacy_mat(i_t) = efficacy; 

fprintf("Percent complete: %d \n",i_t/n_samps*100);   
end

% error_bound_mat
% error_Ahat_mat
% efficacy_mat
xi_rand = xi_rand*2-1;

if nom_opt ~= 1
    save(strcat('LDC_design/rand_',save_label),'error_bound_mat','error_Ahat_mat','efficacy_mat', 'u_lim','nu_lim','xi_rand');
end

end
