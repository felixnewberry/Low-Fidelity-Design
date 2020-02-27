clear all
close all
% clc

save_on = 0; 

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
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These samples suck... 

n_samples = 200; 

load('L_data/Uf_stress')


load('L_data/Uc_stress')
Uc_nom = Uc; 
load('L_data/Uc_stress_opt')
Uc_opt = Uc; 


% load Ub 
load('L_data/Ub_stress')
Ub_nom = U; 
load('L_data/Ub_stress_opt')
Ub_opt = U; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_b_nom = vecnorm(Ub_nom-Uf)./vecnorm(Uf);
error_b_opt = vecnorm(Ub_opt-Uf)./vecnorm(Uf);

n_hist = 20; 
figure
hold on
h1 = histogram(abs(100*error_b_nom),n_hist,'FaceColor',c1);
h2 = histogram(abs(100*error_b_opt),n_hist,'FaceColor',c2);
hold off
legend([h1,h2],{'Nominal','Optimal'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('Relative Error $[\%]$','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf,'Position',size_1)
title('U Mid','Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply transformations + pre-processing. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% try differet r
r = 10; % or 6
n = r+10; 
% r = 6; % or 6
% n = r+2; 

rand_sample = 1:n; 
N = n_samples; 

Uc = Uc_nom; 

% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

1; 

[error_bound_nom,err_Bi_nom, P_s_nom] = my_bound_bi(n,r, A, B, N);

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

% Tansform L 
Phi = A_R* pinv(B_R); 
B = Phi*B; 

B = B/norm(B,'fro'); %? Should I take this step? Probably. 

B_2 = B; 
1; 
% seems to break things... I should run an L shaped sim and save the data
% for further exploration. Plot these errors! The histograms! 

[error_bound_nom2,err_Bi_nom2, P_s_nom2] = my_bound_bi(n,r, A, B, N);

% why is the bound no longer true? What assumptions are being broken? 
% sb after 20 goes to zero. Possibly if things were interpolated this would
% work... 
% Must be breaking some assumption here. What could it be? I've messed with
% B
error_B = norm(A - B)/norm(A); % significant error - about 7 % ... 

% B was previously rank 41, and now it is rank n... This means that sb_n+1
% is always 0... why was this not an issue for the other codes? 

%Brief review of conditions in paper: 

% With this method I interpolate too. Low fidelity goes from 41 to 465... 
% is it better if I do this in advance? I could test this just with nominal
% to see. 

1; 

% 
Uc = Uc_opt; 
% % Normalize matrices for bi-fidelity error bound
B = Uc/norm(Uc,'fro');
[error_bound_opt,err_Bi_opt, P_s_opt] = my_bound_bi(n,r, A, B, N);

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

% Tansform L 
Phi = A_R* pinv(B_R); 
B = Phi*B; 
[error_bound_opt2,err_Bi_opt2, P_s_opt2] = my_bound_bi(n,r, A, B, N);

res_mat =  [error_bound_nom, error_bound_nom2,error_bound_opt,...
    error_bound_opt2; ...
    err_Bi_nom, err_Bi_nom2, err_Bi_opt, err_Bi_opt2]; 

res_tab = array2table(res_mat,'VariableNames',{'Nom','Nom2', 'Opt', 'Opt2'}, 'RowNames',{'Bound','Bi'})

% How is it possible for the bound to break?? 
% Check conditions: 


