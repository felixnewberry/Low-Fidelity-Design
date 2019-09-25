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
%%% Choose optimization method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

check_data = 0; % update with 0 pertubations to bound. 

line_search = 1; 


% I think this is broken for E (=1) just now... 

mode_delta = 1;   
%%% change this
% 0 is nu
% 1 is E
% 2 is corr length
% 3 is sigma
% 4 is theta - changing trapezoid angle of slope of trapezoid. Nominal value 0. Max
% maybe pi/8 to start. (corresponds to 0.4142 in max change to applied load)
% 5 is q: changing mean width of trapezoid

% temporary, try mode_delta = 1 as sigma or corr

% % 2 is h2
% % 3 is h3
% % 4 is h1=h2 and h3

mode_qoi = 0; 
% 0 is line from x = 0, y = 0 to y = 1 change r to 6
% 1 is field 
% 2 is field 


% mode_qoi = 0; 
% mode_qoi = 1; 



% need to fix up with deltas. 
grid_search = 0; % 1515 s when 40x40 grid with 50 repetitions of the bound estimate

random_search = 0; % pce error is about 8 % - stick with grid?
% error improves if I take average of 50 bounds. 
% maybe 5 mintutes? 

plot_tip = 0; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data, vanilla approach first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordinates
load('L_data/dof_coords_c')
dof_coords_c = dof_coords; 
load('L_data/dof_coords_f')
dof_coords_f = dof_coords; 

% %Uf fine
% load('L_data/Uf')
% Uf = U; 
% % Uc course
% load('L_data/Uc')
% Uc = U; 

% xi 
load('fenics_inputs/xi')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% epsilon tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select subset R columns of N samples from A and B. 
% rank of Uc is r=1. 
% section 2.5 says that n should be slightly larger than r. heuristic R = r+10

% Set 
% tol = 1e-4; 

% heuristic n = r+10
r = 10; 
r = 6; % or it may break, not certain on nuances of this. 
% some NaNs for im rank_k_gsqr within matrixID for certain pertubation values. 
% r = 10;
n = r+10; 

% error saturates at about 7 % for r = 10

% r = 5 sees reduction from 4.14% to 3.52 % 
% r = 10 barely changes this. 

% mess with parameters.. 

% nsim = 800; 
nsim = 200; 

n_sample = 200; 
efficacy_vec = zeros(1,n_sample); 

n_bound_reps = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if check_data ==1 

nsim = 200; 


%%% test error bound
X = 0; 

% Find Uh Uc and bound 
% [error_bound,err_Ahat, Uc, Uf] = my_L_low_high(X,nsim, n, r,n_bound_reps, mode_delta, mode_qoi); 


% Find values from field. Also sort, sort is used for current high fidelity
% data. 

[error_bound,err_Ahat,efficacy] = my_L_bound(X,nsim, n, r, mode_delta, n_bound_reps, mode_qoi); 
error_bound
err_Ahat
% Load established Uh and nominal Ul and Ub. 

1; 

if mode_qoi == 0
    load('L_data/Idx_f')
    load('L_data/Idx_c')

    load('L_data/Uf_line')
    Uf = U(Idx_f,:); 
    load('L_data/Uc_line')
    Uc = U(Idx_c,:); 

    load('L_data/x_f')
    load('L_data/x_c')
    load('L_data/Ub_line')
    
    load('L_data/Uc_line_2')
        
        % interploate surface solution
    Uc_int = zeros(length(x_f),nsim); 

    for i_int = 1:nsim
        Uc_int(:,i_int) = interp1(x_c,Uc(:,i_int),x_f); 
        1; 

    end
elseif mode_qoi == 1
    load('L_data/Uf_field')
    Uf = U; 
    load('L_data/Uc_field')
    Uc = U; 
    load('L_data/Ub_field')
    Ub = U; 
    Uc_int = zeros(length(dof_coords_f),nsim); 

    for i_int = 1:nsim
        F = scatteredInterpolant(dof_coords_c(:,1),dof_coords_c(:,2), Uc(:,i_int));
        Uc_int(:,i_int) = F(dof_coords_f(:,1),dof_coords_f(:,2));
    end
    
elseif mode_qoi == 2
    load('L_data/Uf_stress')
    Uf = Uf; 
    load('L_data/Uc_stress')
    Uc = Uc; 
    load('L_data/Ub_stress')
    Ub = U; 
    Uc_int = zeros(length(dof_coords_f),nsim); 

    for i_int = 1:nsim
        F = scatteredInterpolant(dof_coords_c(:,1),dof_coords_c(:,2), Uc(:,i_int));
        Uc_int(:,i_int) = F(dof_coords_f(:,1),dof_coords_f(:,2));
    end
end

if mode_qoi == 0 
    figure
    p1 = plot(x_f,Uf(:,end),'-x','color',c1, 'LineWidth', LW, 'MarkerSize', MS);
    hold on
    p2 = plot(x_c,Uc(:,end),'-o','color',c2, 'LineWidth', LW, 'MarkerSize', MS);
    % p3 = plot(x_f,Uc_int(:,end),'-d','color',c3, 'LineWidth', LW, 'MarkerSize', MS);
    p4 = plot(x_f,Ub(:,end),'-s','color',c3, 'LineWidth', LW, 'MarkerSize', MS);
    hold off
    xlabel('y', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Horizontal displacement', 'interpreter', 'latex', 'fontsize', FS)
    % legend([p1,p2,p3,p4],{'H','L','L_int','B'},'interpreter', 'latex', 'fontsize', FS_leg)
    legend([p1,p2,p4],{'H','L','B'},'interpreter', 'latex', 'fontsize', FS_leg)
    grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
    
    % pdf
    [f_f,x_pdf_f] = ksdensity(Uf(1,:)); 
    [f_c,x_pdf_c] = ksdensity(Uc(1,:)); 
    [f_b,x_pdf_b] = ksdensity(Ub(1,:));

    figure % point 0, 1
    p1 = plot(x_pdf_f,f_f , 'color',c1,'LineWidth',LW);
    hold on
    p2 = plot(x_pdf_c,f_c , 'color',c2,'LineWidth',LW);
    p3 = plot(x_pdf_b,f_b , 'color',c3,'LineWidth',LW);
    hold off
    xlabel('Horizontal displacement, $u(0.0,1.0)$', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('plot of $u(0.0,1.0)$', 'interpreter', 'latex', 'fontsize', FS)
    legend([p1,p2,p3],{'H','L','B'},'interpreter', 'latex', 'fontsize', FS_leg)
    grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
    
    
elseif mode_qoi == 1 || mode_qoi == 2
    figure
    p1 = plot3(dof_coords_f(:,1),dof_coords_f(:,2),Uf(:,end),'x','color',c1);
    hold on
    p2 = plot3(dof_coords_c(:,1),dof_coords_c(:,2),Uc(:,end),'o','color',c2);
    p3 = plot3(dof_coords_f(:,1),dof_coords_f(:,2),Uc_int(:,end),'d','color',c3);
    hold off
    xlabel('x', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('y', 'interpreter', 'latex', 'fontsize', FS)
    zlabel('Horizontal displacement', 'interpreter', 'latex', 'fontsize', FS)
    legend([p1,p2,p3],{'H','L','L_int'},'interpreter', 'latex', 'fontsize', FS_leg)
    grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
    
end


error_L = norm(Uf-Uc_int)/norm(Uf);
error_B = norm(Uf-Ub)/norm(Uf);

fprintf("Low-fidelity error:  %d \n Error Bound:  %d \n Bi-fidelity error:  %d \n",error_L, error_bound, err_Ahat);

% svd
figure
p1 = semilogy(svd(Uf)/max(svd(Uf)),'-x','color',c1,'LineWidth',LW,'MarkerSize',MS); 
hold on
p2 = semilogy(svd(Uc)/max(svd(Uc)),'-o','color',c2,'LineWidth',LW,'MarkerSize',MS); 
% p3 = semilogy(svd(Ub)/max(svd(Ub)),'-x','color',c1,'LineWidth',LW,'MarkerSize',MS); 
xlabel('index $i$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Normalized singular value', 'interpreter', 'latex', 'fontsize', FS)
legend([p1,p2],{'H','L'},'interpreter', 'latex', 'fontsize', FS_leg)
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
xlim([1,30])





end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Line search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if line_search == 1
    
% choose mode to test w, t1, t2, t3. mode [0,1,2,3]
% Update plot according.

% nominal and optimal delta values
% delta_t1_rand = [0, 1.5];
% delta_t3_rand) = [0,29]; 

% delta_vec = -0.5:0.1:0.5;

% nu
% delta_vec = 0.0:0.05:0.4;
% delta_vec = 0.0:0.05:0.5;

% delta_vec = -0.3:0.05:0.3;

% delta_vec = [-0.3,0.3];

% delta_vec = -0.04:0.02:0.04;

% E
% delta_vec = 0:0.025:0.2;
% do finer detail for sigma
delta_vec = -0.5:0.5:0.5;

% corr
% delta_vec = 0.0:0.1:2.6;
% delta_vec = -0.95:0.95:0.95;

% sigma
% delta_vec = -0.05:0.05:0.05;
% delta_vec = -0.95:0.15:0.2;

% delta_vec = -0.95:0.95:0.95;

% Theta
% delta_vec = [-pi/8:pi/16:pi/8];
% delta_vec = [-7*pi/16:pi/16:0];

% delta_q
% delta_vec = [-0.5:0.25:0.5];

if mode_delta == 0
%     delta_vec = -0.5:0.1:0.5;
%     delta_vec = -0.95:0.1:0.5;

    plot_label = '$ \Delta \nu [\%]$';
elseif mode_delta == 1
%     delta_vec = -0.95:0.1:0.5;
    plot_label = '$ \Delta E [\%]$';
elseif mode_delta == 2
    plot_label = '$ \Delta corr length [\%]$';
elseif mode_delta == 3
    plot_label = '$ \Delta \sigma [\%]$';
elseif mode_delta == 4
    plot_label = '$ \Delta \theta [rad] $';    
elseif mode_delta == 5
    plot_label = '$ \Delta q [\%] $';    
end

% delta_vec = -0.5:0.1:0.5;

% t3 
error_bound_mat = zeros(length(delta_vec),1); 

error_Ahat_mat = zeros(length(delta_vec),1);
efficacy_mat = zeros(length(delta_vec),1);

tic

%%% matrixIDvR fails on svd. NaNs in R if maxrnk > rank = 6 of problem... 

for i_test = 1:length(delta_vec)
delta_vec(i_test)
[error_bound,err_Ahat,efficacy] = my_L_bound(delta_vec(i_test),nsim, n, r, mode_delta, n_bound_reps, mode_qoi); 

error_bound_mat(i_test) = error_bound;
error_Ahat_mat(i_test) =  err_Ahat;
efficacy_mat(i_test) = efficacy; 

1; 


end

% if mode == 0
%     save('Beam_design/line_w_1to10','error_bound_mat')
%     save('Beam_design/delta_vec','delta_vec')
% elseif mode == 1
%     save('Beam_design/line_h1_1to10','error_bound_mat')
% elseif mode == 2
%     save('Beam_design/line_h2_1to10','error_bound_mat')
% elseif mode == 3
%     save('Beam_design/line_h3_1to10','error_bound_mat')
% end
% toc

% % Improve contour plot
figure
hold on
p1 = plot(100*delta_vec,100*error_bound_mat,'ob-', 'LineWidth',LW,'MarkerSize',MS); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on

1; 

% Not appropriate to plot Ahat here
figure
hold on
p1 = plot(100*delta_vec,100*error_Ahat_mat,'sr-', 'LineWidth',LW); 
hold off
xlabel(plot_label,'interpreter','latex','Fontsize',FS)
ylabel('Error Bi $[\%]$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on


% load('Beam_design/line_w')
% error_w = error_bound_mat;
% load('Beam_design/line_h1')
% error_h1 = error_bound_mat;
% load('Beam_design/line_h2')
% error_h2 = error_bound_mat;
% load('Beam_design/line_h3')
% error_h3 = error_bound_mat;
% load('Beam_design/delta_vec')
% 
% load('Beam_design/line_w_1to10')
% error_w = error_bound_mat;
% load('Beam_design/line_h1_1to10')
% error_h1 = error_bound_mat;
% load('Beam_design/line_h2_1to10')
% error_h2 = error_bound_mat;
% load('Beam_design/line_h3_1to10')
% error_h3 = error_bound_mat;
% load('Beam_design/delta_vec')

% matlab default colors look good: how to use these? 

% figure
% hold on
% p1 = plot(100*delta_vec,100*error_w,'o-', 'Color',c1,'LineWidth',LW,'MarkerSize',MS); 
% p2 = plot(100*delta_vec,100*error_h1,'x--', 'Color',c2, 'LineWidth',LW,'MarkerSize',MS); 
% p3 = plot(100*delta_vec,100*error_h2,'s-.', 'Color',c3, 'LineWidth',LW,'MarkerSize',MS); 
% p4 = plot(100*delta_vec,100*error_h3,'d:', 'Color',c4, 'LineWidth',LW,'MarkerSize',MS); 
% hold off
% xlabel('$\Delta [\%]$','interpreter','latex','Fontsize',FS)
% ylabel('Error Bound $[\%]$','interpreter','latex','Fontsize',FS)
% legend([p1,p2,p3,p4],{'$w$','$h_1$','$h_2$','$h_3$'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','SouthEast')
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on



end

if grid_search == 1
%%

mode_delta = 4; 

% 40 x 40 takes 196 s
% t3_vec = 5:5:200; 
% t1_vec = [0.005:0.015:0.6];
delta_t3_vec = 0:1:39; 
delta_t1_vec = -0.975:0.075:2; 

% 20 x 20 takes 45 s. 
% t3_vec = 5:10:200; 
% t1_vec = 0.005:0.03:0.6;
% delta_t3_vec = 0:2:39; 
% delta_t1_vec = -0.975:0.15:2; 

tic
error_bound_mat = zeros(length(delta_t1_vec),length(delta_t3_vec)); 
error_Ahat_mat = zeros(length(delta_t1_vec),length(delta_t3_vec)); 
efficacy_mat = zeros(length(delta_t1_vec),length(delta_t3_vec)); 

for i_test = 1:length(delta_t1_vec)
 i_test
for i_t3 = 1:length(delta_t3_vec)

[error_bound,err_Ahat,efficacy] =  my_beam_bound([delta_t1_vec(i_test),delta_t3_vec(i_t3)],nsim, n, r, mode_delta, n_bound_reps);

% results(i_t3,:) = results(i_t3,:) +[error_bound, err_Ahat, efficacy]; 
error_bound_mat(i_test,i_t3) = error_bound; 
error_Ahat_mat(i_test,i_t3) =err_Ahat; 
efficacy_mat(i_test,i_t3) = efficacy; 

% error_efficacy_sum
end

end

toc
% Identify minimum with respect to both A_hat and bound

% A_hat
% Minimum of each column
[min_values, index_t1_A] = min(error_Ahat_mat); 
% Global minimum, column index is t3
[min_value, index_t3_A] = min(min_values); 
index_t1_A = index_t1_A(index_t3_A);
t1_Ahat = delta_t1_vec(index_t1_A);
t3_Ahat = delta_t3_vec(index_t3_A);
min_Ahat = min_value; 

min_bound_Ahat = error_bound_mat(index_t1_A,index_t3_A);

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
min_Ahat_bound = error_Ahat_mat(index_t1_B,index_t3_B);

fprintf('Bound location and values \n');
fprintf('t3: %d, t1=t2: %d \n',t3_bound, t1_bound);
fprintf('Bound: %d \n',min_bound);
fprintf('Ahat: %d \n',min_Ahat_bound);

fprintf('Ahat location and values \n');
fprintf('t3: %d, t1=t2: %d \n',t3_Ahat, t1_Ahat);
fprintf('Bound: %d \n',min_bound_Ahat);
fprintf('Ahat: %d \n',min_Ahat);

figure
hold on
contourf(100*delta_t3_vec,100*delta_t1_vec,100*error_bound_mat)
colorbar
% 5, 0.2
p1 = plot(0,0,'ro','MarkerSize',8,'linewidth',LW);
p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',8,'linewidth',LW);
p3 = plot(100*t3_Ahat,100*t1_Ahat,'rs','MarkerSize',8,'linewidth',LW);
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
contourf(100*delta_t3_vec,100*delta_t1_vec,100*error_Ahat_mat)
colorbar
p1 = plot(0,0,'ro','MarkerSize',8,'linewidth',LW);
p2 = plot(100*t3_bound,100*t1_bound,'rx','MarkerSize',8,'linewidth',LW);
p3 = plot(100*t3_Ahat,100*t1_Ahat,'rs','MarkerSize',8,'linewidth',LW);
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

mode_delta = 4; 
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
error_Ahat_mat = zeros(n_samps,1);
efficacy_mat = zeros(n_samps,1);


% % illustrate best result
delta_t1_rand(1:2) = [0, 1.5];
delta_t3_rand(1:2) = [0,29]; 

for i_t = 1:n_samps
    
[error_bound,err_Ahat,efficacy] =  my_beam_bound([delta_t1_rand(i_t),delta_t3_rand(i_t)],nsim, n, r, mode_delta, n_bound_reps);

% error_bound = ahat_error_est/norm(Uf);
% err_Ahat = norm(Uf-Uf(:,ix)*P_s)/norm(Uf);
% efficacy = error_bound/err_Ahat;

error_bound_mat(i_t) = error_bound;
error_Ahat_mat(i_t) =  err_Ahat;
efficacy_mat(i_t) = efficacy; 

% error_bound
% err_Ahat
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
c_ref_A = psi\error_Ahat_mat(1:N,1); 

% opts = spgSetParms('iterations',10000,'verbosity',0);
opts = spgSetParms('iterations',10000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);

% % solve with spgl1
weights = get_matrix_weights(psi);
Psiweights = psi*weights;
% % sigma is truncation error of PCE (approximated)
sigma =  cross_val_sigma(psi,error_bound_mat(1:N,1));
c_spg = weights*spg_bpdn(Psiweights,error_bound_mat(1:N,1),sigma*norm(error_bound_mat(1:N,1)),opts);

sigma_A =  cross_val_sigma(psi,error_Ahat_mat(1:N,1));
c_spg_A = weights*spg_bpdn(Psiweights,error_Ahat_mat(1:N,1),sigma_A*norm(error_Ahat_mat(1:N,1)),opts);

%%% Start validation

clear psi;

for isim=1:nsim_v
    crow = piset(xi_rand(isim+N,:),index_pc);
    psi(isim,:) = crow(1:P);
end

error_val_ls = norm(error_bound_mat(N+1:end,1)-psi*c_ref)/norm(error_bound_mat(N+1:end,1));
error_val_spg = norm(error_bound_mat(N+1:end,1)-psi*c_spg)/norm(error_bound_mat(N+1:end,1));

error_val_ls_A = norm(error_Ahat_mat(N+1:end,1)-psi*c_ref_A)/norm(error_Ahat_mat(N+1:end,1));
error_val_spg_A = norm(error_Ahat_mat(N+1:end,1)-psi*c_spg_A)/norm(error_Ahat_mat(N+1:end,1));

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






% figure
% contourf(efficacy_sum)
% colorbar
% 
% point_height = 1; 
% figure
% surf(t3_vec,t1_vec,error_Ahat_sum)
% colorbar
% view(2)
% shading interp
% axis tight;
% hold on
% p1 = plot(5,0.2,point_height,'ro','MarkerSize',8,'linewidth',LW);
% p2 = plot(t3_bound,t1_bound,point_height,'rx','MarkerSize',8,'linewidth',LW);
% p3 = plot(t3_Ahat,t1_Ahat,point_height,'rs','MarkerSize',8,'linewidth',LW);
% hold off
% legend([p1,p2,p3],{'Nominal','Bound', '$\hat{A}$ Min'},'interpreter', 'latex', 'fontsize', FS_leg)
% xlabel('t3','interpreter','latex','Fontsize',FS)
% ylabel('t1=t2','interpreter','latex','Fontsize',FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
% title('$\hat{A}$ error','Interpreter','latex')

% figure
% contourf(t3_vec,t1_vec,error_Ahat_sum)
% shading interp
% colorbar





% plot for alireza: 
% By any chance, for the beam problem, do you have the histogram of... 
% the tip displacement errors for the regular and optimized bi-fidelity...
% errors? You mentioned the error reduced to 0.4% from 3% and I was hoping...
% to show in a report the associated error histogram. 

% t3 is 195 
% t2=21 is 0.59
% 
% nominal t3 = 5, t2 = t1 = 0.2

