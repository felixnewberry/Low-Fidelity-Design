%% Test analytic LDC problem 

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
%%% Plot Nominal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find RE... 

load 'x_64.mat'

x_64 = x_64(:,1); 
y_64 = x_64; 

[X,Y] = meshgrid(x_64, y_64); 

u = zeros(size(X)); 
v = zeros(size(X)); 

for i_X = 1:length(X(:))
        u(i_X) = u_ldc(1.0, X(i_X),Y(i_X));
        v(i_X) = v_ldc(1.0, X(i_X),Y(i_X));
end

u_mag = sqrt(u.^2 + v.^2); 

step_quiv = 5; 

figure
[c,h]=contourf(X,Y,u_mag);
set(h, 'edgecolor','none');
hold on

quiver(x_64(1:step_quiv:end),y_64(1:step_quiv:end),u(1:step_quiv:end,1:step_quiv:end), v(1:step_quiv:end,1:step_quiv:end),'Color',c3)

% % plot qoi
% p1 = plot(x_1, y_2,'-.','color',c2,'LineWidth',LW+1);
% p2 = plot(x_1, y_1,'--','color',c3,'LineWidth',LW+1);

% X_arrow = [0.35 0.7];
% Y_arrow = [0.95   0.95];
% hh = annotation('arrow',X_arrow,Y_arrow,'Color','r');
% set(hh, 'LineWidth', LW)

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
% legend([p1,p2],{'U Mid','P Base'},'interpreter', 'latex', 'fontsize', FS_leg)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
set(gcf, 'Position', size_square)
pbaspect([1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate analytic ensemble 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I could calculate a generic Re from nu and u and feed this into the
% analytic model. 

% This model doesn't have pressure singularities. 

%%% Nominal 
Re = 100; 

x_vec = x_64; 
y_p5 = 0.5; 

% v_p5 = zeros(size(x_vec)); 
% u_top = zeros(size(x_vec)); 

% for i_x = 1:length(x_vec)
%     v_p5(i_x) = v_ldc(Re, x_vec(i_x),y_p5);
%     u_top(i_x) = u_ldc(Re, x_vec(i_x),1.0);
% end

v_p5 = v_ldc(Re, x_vec, y_p5); 
u_top = u_ldc(Re, x_vec,1.0);

figure 
plot(x_vec, v_p5)

figure 
plot(x_vec, u_top)

%%% Ensemble 
load u_nu_vec_2.mat
% nu_vec and u_top_vec

Re_vec = u_top_vec./nu_vec;


L_nom = v_ldc(Re_vec, x_vec, y_p5); 

% How to generate pressure? 
% Set pressure in bottom left to be 0. Then implement integration scheme...
% % This would increase the cost alot.


figure
plot(L_nom')

% save('LDC_design/L_nom','L_nom')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analytic Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = u_ldc(Re,x,y)
u = Re*8*((x.^4-2*x.^3+x.^2)*(4*y.^3-2.*y))'; 
end

function v = v_ldc(Re,x,y)
v = Re.*-8*((4*x.^3-6*x.^2+2*x)*(y.^4-y.^2))'; 
end

%%% Re = 1

% function u = u_ldc(x,y)
% u = 8*(x^4-2*x^3+x^2)*(4*y^3-2*y); 
% end
% 
% function v = v_ldc(x,y)
% v = -8*(4*x^3-6*x^2+2*x)*(y^4-y^2); 
% end

