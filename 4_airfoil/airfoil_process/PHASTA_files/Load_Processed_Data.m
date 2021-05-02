function [ NACA_params, NACA_data, NACA_stats ] = ...
    Load_Processed_Data( identifier_string, nfiles, alpha_threshold )

    fprintf('Reading %i files from %s\n',nfiles,identifier_string);
    
    main_data_dir = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
                    'ASEN_6519_Uncertainty_Quantification/Final_Project/data/'];

    if strcmp(identifier_string, 'original-LF') || strcmp(identifier_string, 'original-HF')
        starting_run_number = 1000;
        params_file = [main_data_dir, 'naca_params_posalpha.csv'];
        if strcmp(identifier_string, 'original-LF')
            path = [main_data_dir, 'original-lf/'];
        elseif strcmp(identifier_string, 'original-HF')
            path = [main_data_dir, 'original-hf/'];
        end
    elseif strcmp(identifier_string, 'attack-LF') || strcmp(identifier_string, 'attack-HF')
        starting_run_number = 0;
        params_file = [main_data_dir, 'naca_params_4412alphaneg2to2.csv'];
        if strcmp(identifier_string, 'attack-LF')
            path = [main_data_dir, 'attack-lf/'];
        elseif strcmp(identifier_string, 'attack-HF')
            path = [main_data_dir, 'attack-hf/'];
        end
    elseif strcmp(identifier_string, 'geometry-LF') || strcmp(identifier_string, 'geometry-HF')
        starting_run_number = 0;
        params_file = [main_data_dir, 'naca_params_4412at5percentgeomerr.csv'];
        if strcmp(identifier_string, 'geometry-LF')
            path = [main_data_dir, 'geometry-lf/'];
        elseif strcmp(identifier_string, 'geometry-HF')
            path = [main_data_dir, 'geometry-hf/'];
        end
    elseif strcmp(identifier_string, 'geometryattack-LF') || strcmp(identifier_string, 'geometryattack-HF')
        starting_run_number = 0;
        params_file = [main_data_dir, 'naca_params_4412geomerr10alpha2tom2.csv'];
        if strcmp(identifier_string, 'geometryattack-LF')
            path = [main_data_dir, 'geometryattack-lf/'];
        elseif strcmp(identifier_string, 'geometryattack-HF')
            path = [main_data_dir, 'geometryattack-hf/'];
        end
    elseif strcmp(identifier_string, 'CW-166k-geomerr15pctalpha13p87')
        starting_run_number = 0;
        params_file = [main_data_dir, 'naca_params_4412geomerr15pctalpha13p87.csv'];
        path = [main_data_dir, 'A0a_2face_CW_166k-geomerr15pctalpha13p87/'];
    elseif strcmp(identifier_string, 'CW-166k-geomerr25p20p20palpha13p87')
        starting_run_number = 0;
        params_file = [main_data_dir, 'naca_params_4412geomerr25p20p20palpha13p87.csv'];
        path = [main_data_dir, 'A0a_2face_CW_166k-geomerr25p20p20palpha13p87/'];
    elseif strcmp(identifier_string, 'params_a_170215')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_RANS_204k/16-1-Chef/RunB_4412-AoA-Sweep/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'params_b_170223')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_WMRANS_199k/16-1-Chef/RunA_4412-AoA-Sweep/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'params_b_170223_runs170301')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_WMRANS_174k/16-1-Chef/RunA_4412-AoA-Sweep/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'params_b_170223_runs170330')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_WMRANS_106k/16-1-Chef/RunA_4412-AoA-Sweep/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'params_b_170223_runs170331')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_WMRANS_131k/16-1-Chef/RunA_4412-AoA-Sweep/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'params_c_170411')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_RANS_coarse5x_a/16-1-Chef/RunB_params170411/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_c')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunB_params170411/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunC_params170505/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'RANS_coarse5xA_params_c')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_RANS_coarse5x_a/16-1-Chef/RunB_params170411/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'RANS_coarse5xA_params_d')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_RANS_coarse5x_a/16-1-Chef/RunC_params_170505/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata_convergedMoreMore'];
    elseif strcmp(identifier_string, 'RANS34k_params_d')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_RANS_34k/4-1-Chef/RunA_params_170505/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'RANS14k_params_d')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_RANS_14k/4-1-Chef/RunA_params170505/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'RANS9k_params_d')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_RANS_9k/2-1-Chef/RunA_params170505/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_FUN3D')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/1-1-Chef/RunA_params170505_fun3d/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H01')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH01_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H02')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH02_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H03')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH03_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H04')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH04_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H05')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH05_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H06')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH06_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H07')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH07_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H08')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH08_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H09')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH09_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H10')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH10_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H11')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH11_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H12')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH12_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H13')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH13_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H14')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH14_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    elseif strcmp(identifier_string, 'Euler59k_params_d_H15')
        starting_run_number = 0;
        tmp = '~/NACA_BiFidelity_2D/runs/A0_ntmr_sloped_Euler_59k/8-1-Chef/RunH15_params170505_mappedAoA/';
        params_file = [tmp, 'setting_param_file.dat'];
        path = [tmp, 'pvdata'];
    else
        error('String "%s" not recognized',identifier_string);
    end

    NACA_params = csvread(params_file,1,0);
    NACA_data = Load_NACA_Directory(path, 'wing', nfiles);
    [NACA_params, NACA_data] = ...
        Sort_NACA_Data(NACA_params, NACA_data, starting_run_number, alpha_threshold);
    nfiles = length(NACA_data);
	if alpha_threshold <= 0
        fprintf('No files removed due to alpha constraint, because threshold <= 0\n');
    else
        fprintf('After removing alpha > %.2f, %i files remain\n',alpha_threshold,nfiles);
    end
    
    if nfiles <= 0
        nfiles = length(NACA_data);
    end
    
    npoints = length(NACA_data(1).xnrot);
    
    xnrot_avg  = zeros(npoints,1);
    cp_max = NACA_data(1).cp;
    cp_min = cp_max;
    cp_avg = zeros(npoints,1);
    cp_var = zeros(npoints,1);
    for i = 1:nfiles
        xnrot_avg = xnrot_avg + NACA_data(i).xnrot / nfiles;
        cp_avg = cp_avg + NACA_data(i).cp / nfiles;
        cp_var = cp_var + (NACA_data(i).cp).^2 / nfiles;
        ind = NACA_data(i).cp > cp_max;
        cp_max(ind) = NACA_data(i).cp(ind);
        ind = NACA_data(i).cp < cp_min;
        cp_min(ind) = NACA_data(i).cp(ind);
    end
    cp_var = cp_var - cp_avg.^2;
    cp_std = sqrt(cp_var);
    
    % Create x_avg_twosurf coordinate, which has x<0 on upper surface, x>0 on lower
    % surface, and x=0 on the trailing edge.
    trailing_x = find(xnrot_avg == max(xnrot_avg),1);
    xnrot_avg_twosurf = [xnrot_avg(1:trailing_x); max(xnrot_avg)*2 - xnrot_avg(trailing_x+1:end)];
    xnrot_avg_twosurf = xnrot_avg_twosurf - max(xnrot_avg);
    
    NACA_stats.xnrot_avg = xnrot_avg;
    NACA_stats.cp_max = cp_max;
    NACA_stats.cp_min = cp_min;
    NACA_stats.cp_avg = cp_avg;
    NACA_stats.cp_var = cp_var;
    NACA_stats.cp_std = cp_std;
    NACA_stats.xnrot_avg_twosurf = xnrot_avg_twosurf;

end

