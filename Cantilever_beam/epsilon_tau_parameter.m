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

% type of search 
line_search = 0; 
grid_search = 0; % 1515 s when 40x40 grid with 50 repetitions of the bound estimate
random_search = 0; % pce error is about 8 % - stick with grid?
plot_tip = 1; 

% choose parameters to vary
mode = 5; 

% 0 is w
% 1 is h1
% 2 is h2
% 3 is h3
% 4 is h1 = h2
% 5 is h1=h2 and h3

% error improves if I take average of 50 bounds. 
% maybe 5 minutes? - what If I don't take average... 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data, vanilla approach first
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
% section 2.5 says that n should be slightly larger than r. heuristic R = r+10

% Set 
% tol = 1e-4; 

% heuristic n = r+10
r = 1; 
% r = 10;
% n = r+10; 
n = r+2; 


% nsim = 100; 
nsim = 100; 

n_sample = 100; 
efficacy_vec = zeros(1,n_sample); 

% n_bound_reps = 50; % draw different sets of n samples to compute bound from N low-fidelity realizations
n_bound_reps = 1; % does not make sense to do repetitions - this would require access to more high-fidelity samples. 

% Some questions for Alireza about this. 

% how best to sample? Fix 1:10 and stick with this (only captures some
% realizatios so not the best measure)
% randomly sample 50 times for every point. 
% randomly sample the same 50 samples for each point. - likely this... 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Line search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if line_search == 1
    
% choose mode to test w, t1, t2, t3. mode [0,1,2,3]

% nominal and optimal delta values
% delta_t1_rand = [0, 1.5];
% delta_t3_rand = [0,29]; 

% delta_vec = -0.95:0.1:0.95; 
% delta_vec = -1:0.075:2;
delta_vec = 0:1:40; 


if mode == 0
    plot_label = '$ \Delta w [\%]$';
elseif mode == 1
    plot_label = '$ \Delta h_1 [\%]$';
elseif mode == 2
    plot_label = '$ \Delta h_2 [\%]$';
elseif mode == 3
    plot_label = '$ \Delta h_3 [\%]$';
elseif mode == 4
    plot_label = '$ \Delta h_1 = \Delta h_2 [\%]$';
end

error_bound_mat = zeros(length(delta_vec),1); 
error_Bi_mat = zeros(length(delta_vec),1);
efficacy_mat = zeros(length(delta_vec),1);


for i_test = 1:length(delta_vec)
    
[error_bound,err_Bi,efficacy] =  my_beam_bound(delta_vec(i_test),nsim, n, r, mode, n_bound_reps);

error_bound_mat(i_test) = error_bound;
error_Bi_mat(i_test) =  err_Bi;
efficacy_mat(i_test) = efficacy; 

end


% save line for using 1:N, N = r+1 = 2. 
% if mode == 0
%     save('Beam_design/line_w_pm95','error_bound_mat')
% %     save('Beam_design/delta_vec','delta_vec')
% elseif mode == 1
%     save('Beam_design/line_h1_pm95','error_bound_mat')
% elseif mode == 2
%     save('Beam_design/line_h2_pm95','error_bound_mat')
% elseif mode == 3
%     save('Beam_design/line_h3_pm95','error_bound_mat')
% end

% save('Beam_design/line_h1h2_p2m1','error_bound_mat')
% save('Beam_design/delta_vec_h1h2_p2m1','delta_vec')

% save('Beam_design/line_h3_p35','error_bound_mat')
% save('Beam_design/delta_vec_h3_p35','delta_vec')



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

% 40 x 40 takes 196 s
% t3_vec = 5:5:200; 
% t1_vec = [0.005:0.015:0.6];
delta_t3_vec = 0:1:40; 
delta_t1_vec = -1:0.075:2; 

% 20 x 20 takes 45 s. 
% t3_vec = 5:10:200; 
% t1_vec = 0.005:0.03:0.6;
% delta_t3_vec = 0:2:39; 
% delta_t1_vec = -0.975:0.15:2; 

tic
error_bound_mat = zeros(length(delta_t1_vec),length(delta_t3_vec)); 
error_Bi_mat = zeros(length(delta_t1_vec),length(delta_t3_vec)); 
efficacy_mat = zeros(length(delta_t1_vec),length(delta_t3_vec)); 

for i_test = 1:length(delta_t1_vec)
 i_test
for i_t3 = 1:length(delta_t3_vec)

[error_bound,err_Bi,efficacy] =  my_beam_bound([delta_t1_vec(i_test),delta_t3_vec(i_t3)],nsim, n, r, mode, n_bound_reps);

% results(i_t3,:) = results(i_t3,:) +[error_bound, err_Bi, efficacy]; 
error_bound_mat(i_test,i_t3) = error_bound; 
error_Bi_mat(i_test,i_t3) =err_Bi; 
efficacy_mat(i_test,i_t3) = efficacy; 

% error_efficacy_sum
end

end

toc
% Identify minimum with respect to both A_hat and bound

% A_hat
% Minimum of each column
[min_values, index_t1_A] = min(error_Bi_mat); 
% Global minimum, column index is t3
[min_value, index_t3_A] = min(min_values); 
index_t1_A = index_t1_A(index_t3_A);
t1_Bi = delta_t1_vec(index_t1_A);
t3_Bi = delta_t3_vec(index_t3_A);
min_Bi = min_value; 

min_bound_Bi = error_bound_mat(index_t1_A,index_t3_A);

% Error bound
% Minimum of each column
[min_values, index_t1_B] = min(error_bound_mat); 
% Global minimum, column index is t3
[min_value, index_t3_B] = min(min_values); 
index_t1_B = index_t1_B(index_t3_B);
t1_bound = delta_t1_vec(index_t1_B);
t3_bound = delta_t3_vec(index_t3_B);
min_bound = min_value; 

% min_bound_A = error_est_sum(index_t1_A,index_t3_A); 
min_Bi_bound = error_Bi_mat(index_t1_B,index_t3_B);

fprintf('Bound location and values \n');
fprintf('t3: %d, t1=t2: %d \n',t3_bound, t1_bound);
fprintf('Bound: %d \n',min_bound);
fprintf('Bi: %d \n',min_Bi_bound);

fprintf('Bi location and values \n');
fprintf('t3: %d, t1=t2: %d \n',t3_Bi, t1_Bi);
fprintf('Bound: %d \n',min_bound_Bi);
fprintf('Bi: %d \n',min_Bi);

save('Beam_design/grid_search_N3','delta_t3_vec','delta_t1_vec','error_bound_mat', ...
    'error_Bi_mat','t3_bound', 't1_bound','t3_Bi','t1_Bi','min_bound',...
'min_Bi_bound')

figure
hold on
contourf(100*delta_t3_vec,100*delta_t1_vec,100*error_bound_mat)
colorbar
% 5, 0.2
p1 = plot(0,0,'ro','MarkerSize',8,'linewidth',LW);
p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',8,'linewidth',LW);
p3 = plot(100*t3_Bi,100*t1_Bi,'rs','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Nominal','Optimal', 'Bi'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthWest')
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
p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',8,'linewidth',LW);
p3 = plot(100*t3_Bi,100*t1_Bi,'rs','MarkerSize',8,'linewidth',LW);
hold off
legend([p1,p2,p3],{'Nominal','Optimal', 'Bi'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthWest')
xlabel('$\Delta h_3 [\%]$','interpreter','latex','Fontsize',FS)
ylabel('$\Delta h_1= \Delta h_2$ [\%]','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
title('Bi-fidelity error','Interpreter','latex')


% plot(t3_vec, t1_vec_try,'r-','LineWidth',2)

end

if random_search == 1
%%    

mode = 5; 
n_samps = 200; 

% Identify limits from t1_vec and t3_vec
% t3_vec = 5:5:200; 
% t1_vec = [0.005:0.015:0.6];

delta_t1_vec = [-0.975,2]; 
delta_t3_vec = [0,39];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE fit to random points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N_total = n_samps; 
N = 150; 
nsim_v = N_total - N; 
% Desired polynomial order of PCE
d = 2; 
p = 7; 

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
    
[error_bound,err_Bi,efficacy] =  my_beam_bound([delta_t1_rand(i_t),delta_t3_rand(i_t)],nsim, n, r, mode, n_bound_reps);

% error_bound = Bi_error_est/norm(Uf);
% err_Bi = norm(Uf-Uf(:,ix)*P_s)/norm(Uf);
% efficacy = error_bound/err_Bi;

error_bound_mat(i_t) = error_bound;
error_Bi_mat(i_t) =  err_Bi;
efficacy_mat(i_t) = efficacy; 

% error_bound
% err_Bi
% 1;



end

xi_rand = xi_rand*2-1;

clear psi
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

%%% Start validation

clear psi;

for isim=1:nsim_v
    crow = piset(xi_rand(isim+N,:),index_pc);
    psi(isim,:) = crow(1:P);
end

error_val_ls = norm(error_bound_mat(N+1:end,1)-psi*c_ref)/norm(error_bound_mat(N+1:end,1));
error_val_spg = norm(error_bound_mat(N+1:end,1)-psi*c_spg)/norm(error_bound_mat(N+1:end,1));

error_val_ls_A = norm(error_Bi_mat(N+1:end,1)-psi*c_ref_A)/norm(error_Bi_mat(N+1:end,1));
error_val_spg_A = norm(error_Bi_mat(N+1:end,1)-psi*c_spg_A)/norm(error_Bi_mat(N+1:end,1));

error_ls(i_rep) = error_val_ls; 
error_spg(i_rep) = error_val_spg; 

error_ls_A(i_rep) = error_val_ls_A; 
error_spg_A(i_rep) = error_val_spg_A; 
end

toc
 

% % PCE stats: 
fprintf('PCE Statistics: \n');
fprintf('LS Mean Bound: %d, LS Sd: %d \n',mean(error_ls), std(error_ls));
fprintf('SPG Mean Bound: %d SPG Sd: %d \n',mean(error_spg), std(error_spg));

fprintf('LS Mean Bi: %d, LS Sd: %d \n',mean(error_ls_A), std(error_ls_A));
fprintf('SPG Mean Bi: %d SPG Sd: %d \n',mean(error_spg_A), std(error_spg_A));

% using least squares and N = 150 vs 50 for validation
% p = 7
% 6.9% 
% p = 8, 7.35 %  

% Surface from PCE
%     delta_t3_vec = 5:5:200; 
%     delta_t1_vec = [0.005:0.015:0.6];


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

% nominal t3 = 5, t2 = t1 = 0.2

if plot_tip == 1

mode = 5; 

n_samps = 2; 

results = zeros(n_samps,3); 
error_bound_mat = zeros(n_samps,1); 
error_Bi_mat = zeros(n_samps,1);
efficacy_mat = zeros(n_samps,1);

delta_t1_rand = [0,1.025]; 
delta_t3_rand = [0,29]; 

% tic

% have to pause at end of my_beam_bound 
n_samps = 2; 

for i_t = 1:n_samps
    
[error_bound,err_Bi,efficacy] =  my_beam_bound([delta_t1_rand(i_t),delta_t3_rand(i_t)],nsim, n, r, mode,n_bound_reps);

% error_bound = Bi_error_est/norm(Uf);
% err_Bi = norm(Uf-Uf(:,ix)*P_s)/norm(Uf);
% efficacy = error_bound/err_Bi;

error_bound_mat(i_t) = error_bound; 
error_Bi_mat(i_t) =  err_Bi; 
efficacy_mat(i_t) = efficacy; 

end




end
