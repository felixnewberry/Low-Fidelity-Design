clear all
close all
clc


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
%%% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load coordinates
load('L_data/dof_coords_f')
x_f = dof_coords; 
load('L_data/dof_coords_c')
x_c = dof_coords; 

% load  Uc
load('L_data/Uc_stress')
Uc_nom = Uc; 
load('L_data/Uc_stress_opt')
Uc_opt = Uc; 

% load Ub 
load('L_data/Ub_stress')
Ub_nom = U; 
load('L_data/Ub_stress_opt')
Ub_opt = U; 

% load Uf
load('L_data/Uf_stress')


% Interpolate low-fidelity to high-fidelity mesh
nsim = 200; 
Uc_nom_int = zeros(length(x_f),nsim); 
Uc_opt_int = zeros(length(x_f),nsim); 

for i_int = 1:nsim
    F_nom = scatteredInterpolant(x_c(:,1),x_c(:,2), Uc_nom(:,i_int));
    Uc_nom_int(:,i_int) = F_nom(x_f(:,1),x_f(:,2));
    F_opt = scatteredInterpolant(x_c(:,1),x_c(:,2), Uc_opt(:,i_int));
    Uc_opt_int(:,i_int) = F_opt(x_f(:,1),x_f(:,2));   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relative Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Entire domain
e_L_nom = norm(Uf - Uc_nom_int)/norm(Uf);
e_L_opt = norm(Uf - Uc_opt_int)/norm(Uf);

e_B_nom = norm(Uf - Ub_nom)/norm(Uf);
e_B_opt = norm(Uf - Ub_opt)/norm(Uf);

% local to each point (p)
e_L_nom_p = abs((Uf - Uc_nom_int)./(Uf));
e_L_opt_p = abs((Uf - Uc_opt_int)/(Uf));

e_B_nom_p = abs((Uf - Ub_nom)./(Uf));
e_B_opt_p = abs((Uf - Ub_opt)./(Uf));

%%% what coordinate sees the biggest improvement? 
e_L_nom_points = vecnorm(Uf - Uc_nom_int,2,2)./vecnorm(Uf,2,2);
e_B_nom_points = vecnorm(Uf - Ub_nom,2,2)./vecnorm(Uf,2,2);

e_L_opt_points = vecnorm(Uf - Uc_opt_int,2,2)./vecnorm(Uf,2,2);
e_B_opt_points = vecnorm(Uf - Ub_opt,2,2)./vecnorm(Uf,2,2);

[~,i_points]= max(e_B_nom_points - e_B_opt_points);
 
% from e_B_nom_points(i_points) to e_B_opt_points(i_points)
% what sample sees the biggest improvement?  (16.5 to 4.4)... 

e_L_nom_samples = vecnorm(Uf - Uc_nom_int)./vecnorm(Uf);
e_B_nom_samples = vecnorm(Uf - Ub_nom)./vecnorm(Uf);

e_L_opt_samples = vecnorm(Uf - Uc_opt_int)./vecnorm(Uf);
e_B_opt_samples = vecnorm(Uf - Ub_opt)./vecnorm(Uf);

[~,i_sample]= max(e_B_nom_samples - e_B_opt_samples);
% from e_B_nom_samples(i_sample) to e_B_opt_samples(i_sample) (24.6 to
% 4.18)

% need to be careful I don't select a realization that was used to lift - didn't  
e_B_opt_samples(i_sample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pdf 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % pdf of 200 samples at one point... 
% i_points = 1; 
[f_f,x_pdf_f] = ksdensity(Uf(i_points,:)); 
[f_c_nom,x_pdf_c_nom] = ksdensity(Uc_nom_int(i_points,:)); 
[f_b_nom,x_pdf_b_nom] = ksdensity(Ub_nom(i_points,:));
[f_c_opt,x_pdf_c_opt] = ksdensity(Uc_opt_int(i_points,:)); 
[f_b_opt,x_pdf_b_opt] = ksdensity(Ub_opt(i_points,:));

figure 
p1 = plot(x_pdf_f,f_f , 'color',c1,'LineWidth',LW);
hold on
p2 = plot(x_pdf_c_nom,f_c_nom , 'color',c2,'LineWidth',LW);
p3 = plot(x_pdf_b_nom,f_b_nom , 'color',c3,'LineWidth',LW);
p4 = plot(x_pdf_c_opt,f_c_opt , 'color',c4,'LineWidth',LW);
p5 = plot(x_pdf_b_opt,f_b_opt , 'color',c5,'LineWidth',LW);
hold off
xlabel('Von Mises Stress at point $(0.0, 0.1903)$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('pdf', 'interpreter', 'latex', 'fontsize', FS)
legend([p1,p2,p3,p4,p5],{'H','L Nom','B Nom','L Opt','B Opt'},'interpreter', 'latex', 'fontsize', FS_leg)
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;


% or pdf of all stress for one sample? 
% i_sample = 1; 
[f_f,x_pdf_f] = ksdensity(Uf(:,i_sample)); 
[f_c_nom,x_pdf_c_nom] = ksdensity(Uc_nom_int(:,i_sample)); 
[f_b_nom,x_pdf_b_nom] = ksdensity(Ub_nom(:,i_sample));
[f_c_opt,x_pdf_c_opt] = ksdensity(Uc_opt_int(:,i_sample)); 
[f_b_opt,x_pdf_b_opt] = ksdensity(Ub_opt(:,i_sample));

figure 
p1 = plot(x_pdf_f,f_f , 'color',c1,'LineWidth',LW);
hold on
p2 = plot(x_pdf_c_nom,f_c_nom , 'color',c2,'LineWidth',LW);
p3 = plot(x_pdf_b_nom,f_b_nom , 'color',c3,'LineWidth',LW);
p4 = plot(x_pdf_c_opt,f_c_opt , 'color',c4,'LineWidth',LW);
p5 = plot(x_pdf_b_opt,f_b_opt , 'color',c5,'LineWidth',LW);
hold off
xlabel('Von Mises Stress for realization 30', 'interpreter', 'latex', 'fontsize', FS)
ylabel('pdf', 'interpreter', 'latex', 'fontsize', FS)
legend([p1,p2,p3,p4,p5],{'H','L Nom','B Nom','L Opt','B Opt'},'interpreter', 'latex', 'fontsize', FS_leg)
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Possibly I have to do this in paraview? But possibly not. I could inserta
% white box over that corner... 

% can plot actuall stress, or plot the relative errors. 
% Use Uf(:,30), x_f

figure
plot3(x_f(:,1),x_f(:,2),Uf(:,30),'x')

% Interpolate solution - nah. 
[Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);
Vq = griddata(x_f(:,1),x_f(:,2),Uf(:,30),Xq,Yq); 
figure
surf(Xq, Yq, Vq)
view(0,90)

% Interpolate error nom
[Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);
Vq = griddata(x_f(:,1),x_f(:,2),e_B_nom_p(:,30),Xq,Yq); 
zmax = max(Vq(:)); 

figure
subplot(1,2,1)
surf(Xq, Yq, Vq)
view(0,90)
xlabel('x', 'interpreter', 'latex', 'fontsize', FS)
ylabel('y', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight; 
caxis([0,zmax]); % could do 0.4
% zlim([0,zmax])

% Interpolate error bi
[Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);
Vq = griddata(x_f(:,1),x_f(:,2),e_B_opt_p(:,30),Xq,Yq); 

subplot(1,2,2)
surf(Xq, Yq, Vq)
view(0,90)
xlabel('x', 'interpreter', 'latex', 'fontsize', FS)
ylabel('y', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;
caxis([0,zmax]);
% zlim([0,zmax])
colorbar
% Make grid then make top right region white. Should be pretty easy. 

% % hmm 
% figure
% p1 = plot3(dof_coords_f(:,1),dof_coords_f(:,2),Uf(:,end),'x','color',c1);
% hold on
% p2 = plot3(dof_coords_c(:,1),dof_coords_c(:,2),Uc(:,end),'o','color',c2);
% p3 = plot3(dof_coords_f(:,1),dof_coords_f(:,2),Uc_int(:,end),'d','color',c3);
% hold off
% xlabel('x', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('y', 'interpreter', 'latex', 'fontsize', FS)
% zlabel('Horizontal displacement', 'interpreter', 'latex', 'fontsize', FS)
% legend([p1,p2,p3],{'H','L','L_int'},'interpreter', 'latex', 'fontsize', FS_leg)
% grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;


