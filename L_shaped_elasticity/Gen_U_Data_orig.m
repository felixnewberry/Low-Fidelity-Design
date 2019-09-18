clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;     % Line width
FS_leg = 16; % Font size legend


size_1 = [0,0,670,515]; 
% size_1 = [0, 0, 500, 350]; 

FS = 28;    % Font size axis
FS_axis = 18; 
LW_axis = 2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file generates all data for elliptic problem - both fine and coarse

% Input: Define the grid levels
coarse_level = 1; % 1-4
fine_level = 4;

% Number of samples
nsim = 800;
randn('state',1)

% Define sigma, correlation length, and dimension of the one-D Gaussian defining youngs
% sigma = .33;
% corr_length = .3;
% d = 7;
% youngs_0 = 0.1;
% num_kl_grid = 8 * d; % Note this should be appropriately checked

% sigma = .25;
sigma = .33;
corr_length = 0.3;
d = 7;
youngs_0 = 0.1; % was 0.1. not sure what to make of this. 
num_kl_grid = 8 * d; % Note this should be appropriately checked

% Extract number of degrees of freedom and nodal coordinates
frid = fopen(['./mesh/mesh_' num2str(coarse_level) '.txt'],'r');
mesh_info = fscanf(frid,'%d',[1,2]);
xy_coarse = fscanf(frid,'%g %g',[2 mesh_info(1)]); xy_coarse = xy_coarse';
coarse_grid = size(xy_coarse,1);

frid = fopen(['./mesh/mesh_' num2str(fine_level) '.txt'],'r');
mesh_info = fscanf(frid,'%d',[1,2]);
xy_fine = fscanf(frid,'%g %g',[2 mesh_info(1)]); xy_fine = xy_fine';
fine_grid = size(xy_fine,1);

1; 

%only do one grid, the coarse grid
if fine_grid==coarse_grid
    
    %1 generate xi data
    xi = randn(nsim,d^2);
    fxi = sprintf('fenics_inputs/xi%d.mat', coarse_level);
    save(fxi,'xi','-mat')
    
    %2 generate youngs data
    youngs_samp_gen_gaussian_cov(youngs_0, sigma, corr_length, d, xy_coarse, xi, num_kl_grid);
    
    % Generate mesh file name
    system(['cp ./mesh/mesh_' num2str(coarse_level) '.xml' ' ./mesh/mesh.xml']);
    
    tic;
    %3 run python
    system('sudo python3.6 lshape.py')
    t_c = toc/nsim;
    
    %4 upload data to matrix
    U=zeros(coarse_grid,nsim);
    % when interpolating
%     U=zeros(fine_grid,nsim);

    for i=1:nsim
        filename = [ 'L_data/solution.' num2str(i-1) '.txt' ];
        fileID = fopen(filename,'r');
        U(:,i)=fscanf(fileID, '%f');
        fclose(fileID);
    end
   
    file = sprintf( 'L_data/U%d.mat',coarse_level);
    save(file,'U','-mat')
    system('rm -f L_data/solution.*');
    
else %do same thing, but on coarse and fine grid
    
    
    %1 generate xi data
    xi = randn(nsim,d^2);
    fxi = sprintf('fenics_inputs/xi%d.mat', fine_level);
    save(fxi,'xi','-mat')
    %save 'xi.mat' xi
    
    %%%%%%%%%%%%COARSE GRID
    %2A generate youngs data
    youngs_samp_gen_gaussian_cov(youngs_0, sigma, corr_length, d, xy_coarse, xi, num_kl_grid)
    
    % Generate mesh file name
    system(['cp ./mesh/mesh_' num2str(coarse_level) '.xml' ' ./mesh/mesh.xml'])

    %3A run python
    tic;
    system('sudo python3.6 lshape.py')
    t_c = toc/nsim;
    1; 
    
    %4A upload data to matrix
    U=zeros(coarse_grid,nsim);    
%     U=zeros(fine_grid,nsim);    

    for i=1:nsim
        filename = [ 'L_data/solution.' num2str(i-1) '.txt' ];
        fileID = fopen(filename,'r');
        U(:,i)=fscanf(fileID, '%f');
        fclose(fileID);
    end
    
    file = sprintf( 'L_data/U%d.mat',coarse_level);
    save(file,'U','-mat')

    %%%%%%%%%%%%FINE GRID
    %2B generate youngs data
    youngs_samp_gen_gaussian_cov(youngs_0, sigma, corr_length, d, xy_fine, xi, num_kl_grid)
    
    % Generate mesh file name
    system(['cp ./mesh/mesh_' num2str(fine_level) '.xml' ' ./mesh/mesh.xml'])

    %3B run python
    tic;
    system('sudo python3.6 lshape.py')
    t_f = toc/nsim;
    
    %4A upload data to matrix    
    U=zeros(fine_grid,nsim);   
    for i=1:nsim
        filename = [ 'L_data/solution.' num2str(i-1) '.txt' ];
        fileID = fopen(filename,'r');
        U(:,i)=fscanf(fileID, '%f');
        fclose(fileID);
    end
    
    file = sprintf( 'L_data/U%d.mat',fine_level);
    save(file,'U','-mat')
    system('rm -f L_data/solution.*');
    
end
% I think the xy_fine values correspond to Youngs but not to the solution
% field U. 

1; 

load('L_data/U4')
load('L_data/dof_coords_f')

dof_coords(450,:)

% histogram
figure % point 0, 1
histogram(U(450,:),'Normalization','probability')
xlabel('Horizontal displacement, $u(0.0,1.0)$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('plot of $u(0.0,1.0)$', 'interpreter', 'latex', 'fontsize', FS)
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;

% pdf
[f,x_pdf] = ksdensity(U(450,:)); 
figure % point 0, 1
plot(x_pdf,f ,'k','LineWidth',LW)
xlabel('Horizontal displacement, $u(0.0,1.0)$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('plot of $u(0.0,1.0)$', 'interpreter', 'latex', 'fontsize', FS)
grid on; set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on; axis tight;


figure
plot3(dof_coords(:,1),dof_coords(:,2),U(:,end),'x')

