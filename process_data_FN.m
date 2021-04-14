clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Retrive CSV data, compute Cp and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[NparL, NdatL, NstatL] = Load_Processed_Data(id_str, nfiles, alpha_threshold);


% Step through 500 runs... 

directory_path = "../../RunA_felix_low/**/*wing.csv";
% directory_path = "../../RunA_felix_high/**/*-convergedMoreMore*wing.csv";
% hypen or -? Likely doesn't matter. 
% List all files.
files = dir(directory_path);

1; 


% This is effective. Yields name that reads either as 
% pvdata_params-c170505-499_wing.csv 
% or 
% pvdata_params-c170505-499_probe.csv 

% Just need to use a general load and it should be okay. 

% Hopefully this leaves files as the right size... 

% Load_NACA_Data calls Load_Paraview_Data ... 


% Want something like "../../RunA_felix_low/Run-001/pvdata_params-c-170505-001_wing.csv"

1; 

for i = 1:length(files)
    directory_path = files(i).folder;
    directory_path = strcat(directory_path,'/');
    1; 
    
    % needs directory path and data_filename
    [NACA_data(i).IDstr, ...
             NACA_data(i).x, ...
             NACA_data(i).y, ...
             NACA_data(i).p, ...
             NACA_data(i).cp, ...
             NACA_data(i).ref] = Load_NACA_Data(directory_path, files(i).name);
    NACA_data(i).filename = files(i).name;
    % not sure what the .name does... 
    
end


% Does this load everything to NACA_data? Save if it does. 
% I now have relevant data to do something... 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sort data so that it progresses about circumference of wing. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a = uigetfile('*.dat','Select a file')
NACA_params = importdata("setting_param_file.dat",',',1);
NACA_params = NACA_params.data;

% load("setting_param_file.dat")
% NACA_params = csvread(params_file,1,0);


% NACA_data = Load_NACA_Directory(path, 'wing', nfiles);
starting_run_number = 0; 
alpha_threshold = 7; 

1; 

1; 

[NACA_params, NACA_data] = ...
    Sort_NACA_Data(NACA_params, NACA_data, starting_run_number, alpha_threshold);
nfiles = length(NACA_data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot data to check 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pluck out x to plot
% extract Cp and x

1; 

x = zeros(length(NACA_data(1).cp),500); 
Cp = zeros(length(NACA_data(1).cp),500); 


for i = 1:length(files)
    x(:,i) = NACA_data(i).x_twosurf;
    Cp(:,i) = NACA_data(i).cp; 
end

% Cp is just a vector - want it to be a matrix... 
figure
plot(x,Cp)
ax = gca; 
ax.YDir = 'reverse'

save('Cp_low_Cb1_mp2','Cp')
%save('Cp_low_','Cp')
% save('Cp_high','Cp')
% save('x_high','x')
% I have low and high. Now make this a smooth progression to matrix... ie setup and run... 
    
% Also test current low-fidelity efficacy in new script. Streamline process later. % send this to laptop. 
