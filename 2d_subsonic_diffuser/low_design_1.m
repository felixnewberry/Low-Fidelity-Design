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
%%% Nominal data - plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_total = 100; 

load uniform_500.mat;
assembledData_2 = load('assembledDwall_500');
% Use upper and lower area together
Uf = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,4));

load('low_nom.mat')
Uc = sqrt(run_data(:,1).^2+run_data(:,2).^2); 


% Plot 

error_L_vec = abs(Uf - Uc)./abs(Uf);
error_L = norm(Uf-Uc)/norm(Uc); 

figure
p1 = plot(Uf,'o','color',c1,'MarkerSize',MS, 'LineWidth', LW);
hold on
p2 = plot(Uc,'x','color',c2,'MarkerSize',MS, 'LineWidth', LW);
hold off
xlabel('Sample $i$','interpreter','latex','Fontsize',FS)
ylabel('Recirculation length','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1,p2],{'H','Nominal L'},'interpreter', 'latex', 'fontsize', FS_leg)
% title('Mid Velocity Realization','Interpreter','latex')
set(gcf,'Position',size_large)


% I wonder if the crappy samples have some correlation in frequency?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process low-fid AIP QoI - to be moved to simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% New QoI

vL = importdata('v_aip.txt');
uL = importdata('u_aip.txt');
pL = importdata('p_aip.txt');
probe = importdata('p_probe.txt');
prL = pL./probe; 

% Plot data

figure
plot(vL)

figure
plot(uL)

figure
plot(pL)

figure
plot(prL)

% Save low-fid_data - would be better if this was pre-computed in sample. 
save('low_nom','vL','uL','pL','prL','probe')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process high-fid AIP QoI - to be moved to simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Nominal bi-fidelity 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% r = 5; 
% n = r+10; 
% 
% % Subset of vectors for bi-fidelity error estimate
% rand_sample = 1:n; 
% 
% B = Uc/norm(Uc,'fro');
% A = Uf/norm(Uf,'fro');
% 
% B_R = B(:,rand_sample);
% A_R = A(:,rand_sample);
% 
% % Obtain column skeleton of P
% [P_s,ix] = matrixIDvR(B,r);
% 
% % Error bound inputs
% normC = norm(P_s);
% sb = svd(B); 
% err_Bhat = norm(B-B(:,ix)*P_s); 
% N = 100;
% 
% % % % Compute epsilon tau... 
% [~, ahat_error_est,~, ~,~] = ...
%     mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,n);
% 
% error_bound_vec = ahat_error_est/norm(A);
% 
% err_Bi_vec = norm(A-A(:,ix)*P_s)/norm(A);
% err_low_vec = norm(Uc - Uf)/norm(Uf); 


% Sample is a point sample... Hmm. Recirculating region. Hmm. 
% Not amenable to bi-fidelity approach. 

% Instead I could predict the pressure of the line at the intake plane? 
% Check Ryan skinner's work and see if I have this data available.
% I can probably still compute this from the restart... Just dont have time
%-averaged or convergence plots so I would have to recheck this. 

% a bi-fidelity tractable QoI - ie vector or field data. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process AIP QoI - to be moved to simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


