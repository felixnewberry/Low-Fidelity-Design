close all
clear all
clc

% Then apply to both lifting technique and basis reduction and compare
% performance. 
% Are these easy to compare? 

% epsilon tau: max eigenvalue of the difference between H and tau L
% grammians. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;     % Line width
MS = 8;     % Marker Size
FS_leg = 16; % Font size legend

size_1 = [0,0,670,515]; 
size_2 = [0,0,1340,515]; 

% size_1 = [0, 0, 500, 350]; 

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
%%% Choose optimization method and mode 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

point_test = 0; % update with 0 pertubations to bound. 
individual_sensitivity = 0; 

% I think this is broken for E (=1) just now... 

mode_delta = 6; 

%%% change this
% 0 is nu
% 1 is E
% 2 is corr length
% 3 is sigma
% 4 is theta - changing trapezoid angle of slope of trapezoid. Nominal value 0. Max
% maybe pi/8 to start. (corresponds to 0.4142 in max change to applied load)
% 5 is q: changing mean width of trapezoid
% 6 nu, corr, sigma and theta

% need to fix up with deltas. 
% grid_search = 0; % 1515 s when 40x40 grid with 50 repetitions of the bound estimate

random_search = 0; % pce error is about 8 % - stick with grid?
% maybe 5 mintutes? 
nom_opt = 1; 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordinates
load('L_data/dof_coords_c')
dof_coords_c = dof_coords; 
load('L_data/dof_coords_f')
dof_coords_f = dof_coords; 

% xi 
load('fenics_inputs/xi')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% epsilon tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% heuristic n = r+10
r = 6; 
% r = 10;
n = r+10; 

% error saturates at about 7 % for r = 10

% r = 5 sees reduction from 4.14% to 3.52 % 
% r = 10 barely changes this. 

% mess with parameters.. 

nsim = 800; 
% nsim = 200; 

% n_sample = 200; 

n_bound_reps = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if point_test ==1 


%%% test error bound
% Nominal
X = [0; 0; 0; 0]; 
% Optimal stress
% X = [0.3027; 2.3445; -0.4413; 0.3113]; 


[error_bound,error_Bi] = my_L_bound(X,nsim, n, r, mode_delta, nom_opt); 

% Load established Uh and nominal Ul and Ub. 

% save('L_design/point_test_nom','error_bound','error_Bi');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Line search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if individual_sensitivity == 1
    
% mode_delta:
% 0 is nu
% 1 is E
% 2 is corr length
% 3 is sigma
% 4 is theta - changing trapezoid angle of slope of trapezoid. Nominal value 0. Max
% maybe pi/8 to start. (corresponds to 0.4142 in max change to applied load)
% 5 is q: changing mean width of trapezoid

if mode_delta == 0
    %%% 0: nu 
    delta_vec = 0.0:0.1:0.7;
    save_label = 'L_design/individual_qoi_nu_0_p7';
elseif mode_delta == 1
    %%% 1: E
    delta_vec = -0.5:0.25:0.5;
    save_label = 'L_design/individual_qoi_E_mp5_p5'; 
elseif mode_delta == 2  
    %%% 2: Corr
    delta_vec = 0.0:0.2:2.6;
    save_label = 'L_design/individual_qoi_cor_0_2p6';
    % delta_vec = 0.0:0.05:0.2;
%     save_label = 'L_design/individual_qoi_cor_0_p2';
elseif mode_delta == 3
    %%% 3: Sigma
    delta_vec = -1:0.1:0;
    save_label = 'L_design/individual_qoi_sig_mp8_0';
elseif mode_delta == 4
    %%% 4: Theta
%     delta_vec = [-pi/8:pi/32:pi/8];
%     save_label = 'L_design/individual_qoi_theta_m_pi_o8_pi_o8'; 
    delta_vec = [0:pi/32:pi/4];
    save_label = 'L_design/individual_qoi_theta_m_0_pi_o4'; 
elseif mode_delta == 5
    %%% 5: q
    delta_vec = [-0.75:0.25:0.75];
    save_label = 'L_design/individual_qoi_q_mp75_p75';
end

% % % First test, does parameter display bound sensitivity? 
% delta_vec = -0.5:0.1:0.5;
% save_label = strcat('L_design/individual_qoi_',num2str(mode_delta),'_mp5_p5');

% Store all three qoi for all sample points
error_bound_mat = zeros(length(delta_vec),3); 
error_Bi_mat = zeros(length(delta_vec),3);

tic

%%% matrixIDvR fails on svd. NaNs in R if maxrnk > rank = 6 of problem... 

for i_test = 1:length(delta_vec)

[error_bound,error_Bi] = my_L_bound(delta_vec(i_test),nsim, n, r, mode_delta, nom_opt); 

error_bound_mat(i_test,:) = error_bound;
error_Bi_mat(i_test,:) =  error_Bi;
end
 

save(save_label,'error_bound_mat','error_Bi_mat', 'delta_vec')

figure
hold on;
p1 = plot(100*delta_vec,100*error_bound_mat(:,1),'o-', 'Color',c1, 'LineWidth',LW,'MarkerSize',MS); 
p2 = plot(100*delta_vec,100*error_bound_mat(:,2),'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
p3 = plot(100*delta_vec,100*error_bound_mat(:,3),'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
hold off; 
legend([p1,p2,p3],{'$d$ Line','$d$ Field', '$\sigma$ Field'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('$\Delta$','interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on


end

if random_search == 1
%%    

n_samps = 200; 

% For s field
delta_nu_vec = [0,0.7]; 
delta_corr_vec = [1.8,2.6];
delta_sigma_vec = [-0.8,0];
delta_theta_vec = [-pi/8,pi/4];

% % for d line and field
% delta_nu_vec = [0,0.7]; 
% delta_corr_vec = [-0.3,0.6];
% delta_sigma_vec = [0,0];
% delta_theta_vec = [-pi/8,pi/4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE fit to random points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N_total = n_samps; 
N = 150; 
nsim_v = N_total - N; 
% Desired polynomial order of PCE
d = 4; 
p = 3; 

% index_pc = nD_polynomial_array
% assembles index of polynomials for a given stochastic dimension and 
% polynomail order. 

% one column for each variable (ie dimension d). Number of rows for set of
% P basis functions based on order pol (or p) ie 0 0, 1 0, 0 1 for p = 1
% yielding P = 3

index_pc = nD_polynomial_array(d,p); 

P = size(index_pc,1);


psi = zeros(N,P);

n_pc_reps = 1; 
error_ls = zeros(n_pc_reps,1); 
error_spg = zeros(n_pc_reps,1); 

tic
for i_rep = 1:n_pc_reps

    
xi_rand = rand(n_samps,4);

delta_nu_rand = delta_nu_vec(1)+ (delta_nu_vec(2)-delta_nu_vec(1)).*xi_rand(:,1); 
delta_corr_rand = delta_corr_vec(1)+ (delta_corr_vec(2)-delta_corr_vec(1)).*xi_rand(:,2); 
delta_sig_rand = delta_sigma_vec(1)+ (delta_sigma_vec(2)-delta_sigma_vec(1)).*xi_rand(:,3); 
delta_theta_rand = delta_theta_vec(1)+ (delta_theta_vec(2)-delta_theta_vec(1)).*xi_rand(:,4); 

delta_mat = [delta_nu_rand'; delta_corr_rand'; delta_sig_rand'; delta_theta_rand'];

% results = zeros(n_samps,3);

error_bound_mat = zeros(n_samps,3); 
error_Bi_mat = zeros(n_samps,3);


% % % illustrate best result
% delta_nu_rand(1:2) = [0, 1.5];
% delta_corr_rand(1:2) = [0,29]; 

% n_samps = 1; 
% r = 10     
% min pc: 0, 2.6, -0.95 - double check this. 
% min sample: 0.0334, 2.1908, -0.9392
% delta_mat = [0, 0.0, 0,0.0334; 0, 2.6, 2.4, 2.1908; 0, -0.95, 0, -0.9392]; 
% 
% delta_mat = [0.2274; 2.3557; -0.7113]; 
% 
% delta_mat = [0; 0; 0]; 
% r = 10 minimum value is -3.0334 (evidently unrealistic) error of
% 6.9 %. 
% Minimum is corner position at 0, 2.6, -0.95... 
% This seems fairly useless. Check anway... 
% minimum sample value was bound 0.48 with bi 2.13 at [0.0334, 2.1908,
% -0.9392

% r = 6: 
% delta_mat = [0, 0.2274, 0,0.3692; 0, 2.3557, 2.4, 1.8; 0, -0.7308, 0, -0.7308]; 

% min pce: 0.2274, 2.3557, -0.7113
% min sample:  0.3692, 1.8, -0.7308 

for i_t = 1:n_samps
    
[error_bound,error_Bi] = my_L_bound(delta_mat(:,i_t),nsim, n, r, mode_delta, nom_opt); 

error_bound_mat(i_t,:) = error_bound;
error_Bi_mat(i_t,:) =  error_Bi;
i_t
error_bound
error_Bi
1;
end

error_bound_mat(1:n_samps)
error_Bi_mat(1:n_samps)


% save('L_design/nu_corr_sigma_theta_r6','error_bound_mat','error_Bi_mat','delta_mat','xi_rand','delta_nu_vec','delta_corr_vec','delta_sigma_vec');

% save('L_design/nu_corr_sigma_theta_r6','error_bound_mat','error_Bi_mat',...
%     'delta_mat','xi_rand','delta_nu_vec','delta_corr_vec','delta_sigma_vec', 'delta_theta_vec');
% 
% save('L_design/nu_corr_theta_r6_d','error_bound_mat','error_Bi_mat',...
%     'delta_mat','xi_rand','delta_nu_vec','delta_corr_vec','delta_sigma_vec', 'delta_theta_vec');



end

end

if nom_opt == 1
    
    delta_mat = [zeros(1,4); 0.0576, 2.3855, -0.3868, 0.1075]'; 
    

    for i_t = 1:2
        [error_bound,error_Bi] = my_L_bound(delta_mat(:,i_t),nsim, n, r, mode_delta, nom_opt); 
    end

    
end

