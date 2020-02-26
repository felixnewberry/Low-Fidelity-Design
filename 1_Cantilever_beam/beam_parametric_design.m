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
random_search = 0; % pce error is 15 % - not great
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
% delta_test = 1; % +-0.95

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

nsim = 3000; 
% Test nsim of 3000? % 9 21 to 9 23  for 100, 9 24 to 9 54. 
% 3000: 

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

toc

% % Identify minimum with respect to both A_hat and bound
% 
% % A_hat
% % Minimum of each column
% [min_values, index_t1_A] = min(error_Bi_mat); 
% % Global minimum, column index is t3
% [min_value, index_t3_A] = min(min_values); 
% index_t1_A = index_t1_A(index_t3_A);
% t1_Bi = delta_t1_vec(index_t1_A);
% t3_Bi = delta_t3_vec(index_t3_A);
% min_Bi = min_value; 
% 
% min_bound_Bi = error_bound_mat(index_t1_A,index_t3_A);
% 
% % Error bound
% % Minimum of each column
% [min_values, index_t1_B] = min(error_bound_mat); 
% % Global minimum, column index is t3
% [min_value, index_t3_B] = min(min_values); 
% index_t1_B = index_t1_B(index_t3_B);
% t1_bound = delta_t1_vec(index_t1_B);
% t3_bound = delta_t3_vec(index_t3_B);
% min_bound = min_value; 
% 
% % min_bound_A = error_est_sum(index_t1_A,index_t3_A); 
% min_Bi_bound = error_Bi_mat(index_t1_B,index_t3_B);
% 
% fprintf('Bound location and values \n');
% fprintf('t3: %d, t1=t2: %d \n',t3_bound, t1_bound);
% fprintf('Bound: %d \n',min_bound);
% fprintf('Bi: %d \n',min_Bi_bound);
% 
% fprintf('Bi location and values \n');
% fprintf('t3: %d, t1=t2: %d \n',t3_Bi, t1_Bi);
% fprintf('Bound: %d \n',min_bound_Bi);
% fprintf('Bi: %d \n',min_Bi);

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

if random_search == 1
%%    

mode = 5; 
n_samps = 200; 

delta_t1_vec = [-0.975,2]; 
delta_t3_vec = [0,40];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE fit to random points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_total = n_samps; 
N = 200; 
nsim_v = N_total - N; 

tic

xi_rand = rand(n_samps,2);

delta_t1_rand = delta_t1_vec(1)+ (delta_t1_vec(2)-delta_t1_vec(1)).*xi_rand(:,1); 
delta_t3_rand = delta_t3_vec(1)+ (delta_t3_vec(2)-delta_t3_vec(1)).*xi_rand(:,2); 
    
results = zeros(n_samps,3); 
error_bound_mat = zeros(n_samps,1); 
error_Bi_mat = zeros(n_samps,1);
efficacy_mat = zeros(n_samps,1);

% % illustrate best result
delta_t1_rand(1:2) = [0, 1.5];
delta_t3_rand(1:2) = [0,29]; 

for i_t = 1:n_samps
    
[error_bound,err_Bi,efficacy] =  my_beam_bound([delta_t1_rand(i_t),delta_t3_rand(i_t)],nsim, n, r, mode);

error_bound_mat(i_t) = error_bound;
error_Bi_mat(i_t) =  err_Bi;
efficacy_mat(i_t) = efficacy; 
end

% Desired polynomial order of PCE
d = 2; 
p = 7; 

index_pc = nD_polynomial_array(d,p); 
% assembles index of polynomials for a given stochastic dimension and 
% polynomail order. 
% one column for each variable (ie dimension d). Number of rows for set of
% P basis functions based on order pol (or p) ie 0 0, 1 0, 0 1 for p = 1
% yielding P = 3
P = size(index_pc,1);

xi_rand = xi_rand*2-1;

clear psi
psi = zeros(N,P);
tic
for isim=1:N
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_ref = piset(xi_rand(isim,:),index_pc);
    psi(isim,:) = crow_ref(1:P);
end

% Solve with least squares
c_ref = psi\error_bound_mat(1:N,1); 
c_ref_A = psi\error_Bi_mat(1:N,1); 

% opts = spgSetParms('iterations',10000,'verbosity',0);
opts = spgSetParms('iterations',10000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);

% % solve with spgl1
weights = get_matrix_weights(psi);
Psiweights = psi*weights;
% % sigma is truncation error of PCE (approximated)
sigma =  cross_val_sigma(psi,error_bound_mat(1:N,1));
c_spg = weights*spg_bpdn(Psiweights,error_bound_mat(1:N,1),sigma*norm(error_bound_mat(1:N,1)),opts);

sigma_A =  cross_val_sigma(psi,error_Bi_mat(1:N,1));
c_spg_A = weights*spg_bpdn(Psiweights,error_Bi_mat(1:N,1),sigma_A*norm(error_Bi_mat(1:N,1)),opts);

% %%% Start validation % if considering LS
% 
% clear psi;
% psi = zeros(nsim_v,P);
% for isim=1:nsim_v
%     crow = piset(xi_rand(isim+N,:),index_pc);
%     psi(isim,:) = crow(1:P);
% end

% error_ls = norm(error_bound_mat(N+1:end,1)-psi*c_ref)/norm(error_bound_mat(N+1:end,1));
% error_spg = norm(error_bound_mat(N+1:end,1)-psi*c_spg)/norm(error_bound_mat(N+1:end,1));
error_spg = sigma; 
% error_ls_A = norm(error_Bi_mat(N+1:end,1)-psi*c_ref_A)/norm(error_Bi_mat(N+1:end,1));
% error_spg_A = norm(error_Bi_mat(N+1:end,1)-psi*c_spg_A)/norm(error_Bi_mat(N+1:end,1));
error_spg_A = sigma_A;

toc
 
% % PCE stats: 
fprintf('PCE Errors: \n');
% fprintf('LS Bound: %d \n', error_ls);
fprintf('SPG Bound: %d \n', error_spg);

% fprintf('LS Bi: %d \n', error_ls_A);
fprintf('SPG Mean Bi: %d \n', error_spg_A);

% using least squares and N = 150 vs 50 for validation
% p = 7
% 15 % - not great... 


delta_t1_vec = linspace(-0.975,2,40); 
delta_t3_vec = linspace(0,39,40);

[X,Y] = meshgrid(delta_t1_vec,delta_t3_vec);

[xx,yy] = size(X); 
% Transfrom to vector: 
XX = X(:); YY = Y(:); 

% Want -1 to 1 RV. 
%     t1_lim = [0.005,0.6]; 
%     t3_lim = [5,200];

range_t1 = delta_t1_vec(end)-delta_t1_vec(1); 
range_t3 = delta_t3_vec(end)-delta_t3_vec(1); 

XX_xi = (XX-delta_t1_vec(1))/range_t1*2-1;
YY_xi = (YY-delta_t3_vec(1))/range_t3*2-1;
xi_grid = [XX_xi,YY_xi]; 

clear psi 

for isim=1:length(XX)
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_ref = piset(xi_grid(isim,:),index_pc);
    psi(isim,:) = crow_ref(1:P);
end

% estimate response surface
if error_ls <= error_spg
    ZZ = psi*c_ref;
else
    ZZ = psi*c_spg;
end

if error_ls_A <= error_spg_A
    ZZ_A = psi*c_ref_A;
else
    ZZ_A = psi*c_spg_A;
end


% Reshape as matrix
Z = reshape(ZZ,[xx,yy]);
Z_A = reshape(ZZ_A,[xx,yy]);


% interesting. Change range? See if alternative pce 

    % Error bound
% Minimum of each column
[min_values, index_t1_B] = min(Z); 
% Global minimum, column index is t3
[min_value, index_t3_B] = min(min_values); 
index_t1_B = index_t1_B(index_t3_B);
t1_bound = delta_t1_vec(index_t1_B);
t3_bound = delta_t3_vec(index_t3_B);

figure
hold on
contourf(100*Y,100*X,100*Z)
colorbar
% 5, 0.2
p1 = plot(0,0,'ro','MarkerSize',MS,'linewidth',LW);
p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',MS,'linewidth',LW);
hold off
legend([p1,p2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('$\Delta h_3 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1 = \Delta h_2 [\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
title('Bound error','Interpreter','latex')
% title('Error bound','Interpreter','latex')
% to make comparison to earlier contour
%     caxis([0.04,0.35])

figure
hold on
contourf(100*Y,100*X,100*Z_A)
colorbar
% 5, 0.2
p1 = plot(0,0,'ro','MarkerSize',MS,'linewidth',LW);
p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',MS,'linewidth',LW);
hold off
legend([p1,p2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('$\Delta h_3 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1 = \Delta h_2 [\%]$','interpreter','latex','Fontsize',FS)
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

