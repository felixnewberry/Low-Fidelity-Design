function [ naca_params, data_out ] = ...
    Sort_NACA_Data( naca_params, data, starting_run_number, thresh )
%%%
% Sorts data points within each entry of the NACA data structure, such that they progress
% around the circumference of the deformed wing. This is done by computing the camber line
% of the deformed airfoil, and identifying nodes above and below it. These two sets of
% nodes are then sorted by increasing and decreasing x-coordinate, respectively.
%%%

    
    try
        
        load_str = 'Sorting nodal point data...';
        hwait = waitbar(0,sprintf('%s %i / %i',load_str,1,length(data)), ...
                        'CreateCancelBtn','setappdata(gcbf,''cancel_loading'',1)');
        setappdata(hwait,'cancel_loading',0);
        
        n_skipped = 0;
        run_numbers_kept = [];
        for i = 1:length(data)
            
            if getappdata(hwait,'cancel_loading')
                delete(hwait);
                error('Nodal point sorting canceled by user. Program stopped.');
            end
            waitbar(i/length(data),hwait,sprintf('%s %i / %i',load_str,i,length(data)));
            
            if isnan(starting_run_number)
                % Specify NaN starting run number if loading a single run.
                run_number = 1;
            else
                % Extract the run number from the loaded file's ID string (from its filename).
                run_number = strsplit(data(i).IDstr,'-');
                run_number = str2num(run_number{end}); %#ok<ST2NM>
                run_number = run_number - starting_run_number;
            end
            
            1;
            
            % Grab the parameters used to generate this realization of airfoil data.
            m = naca_params(run_number,1);
            p = naca_params(run_number,2);
            t = naca_params(run_number,3);
            a = naca_params(run_number,4);
            c = naca_params(run_number,5);
            
            % Skip if angle of attack is above a certain threshold, otherwise do add run
            % number to list of runs we're keeping.
            if thresh > 0 && a > thresh
                n_skipped = n_skipped + 1;
                continue;
            else
                run_numbers_kept = [run_numbers_kept, run_number];
            end
            
            % Rotate surface coordinate pairs back to 0 angle-of-attack.
            arad = a * pi / 180;
            rotmat = [cos(arad), -sin(arad); sin(arad), cos(arad)];
            tmp = rotmat * [data(i).x,data(i).y]';
            xrot = tmp(1,:)';
            yrot = tmp(2,:)';
            
            % Normalize rotated x-values.
            xnrot = xrot - min(xrot);
            xnrot = xnrot / max(xnrot);
            
            % Sort x-coordinates.
            [xrot_sorted, order] = sort(xrot);
            
            % Calculate the camber line for this airfoil, rotated to zero AoA.
            surface_sorting_line = NACA_Camber_Line(xrot_sorted, m, p, 0.0, c);
            
            % Camber line doesn't work for points where x < 0, so fix for this case with
            % some wiggle room.
            tmp_indices = xrot_sorted < 0.1;
            [~,ind] = min(xrot_sorted);
            surface_sorting_line(tmp_indices) = yrot(order(ind));
            
            % Determine what points are on top or bottom surface of airfoil.
            top_indices = yrot(order) > surface_sorting_line;
            bot_indices = ~top_indices;
            
            % Split the order so top indices go first, followed by negative indices in
            % reverse.
            order = [order(top_indices); flip(order(bot_indices))];
            
            % Re-order the elements of the data matrix corresponding to this run.
            data_out(i-n_skipped).IDstr = data(i).IDstr; %#ok<*AGROW>
            data_out(i-n_skipped).x     = data(i).x(order);
            data_out(i-n_skipped).xrot  = xrot(order);
            data_out(i-n_skipped).xnrot = xnrot(order);
            data_out(i-n_skipped).y     = data(i).y(order);
            data_out(i-n_skipped).yrot  = yrot(order);
            data_out(i-n_skipped).p     = data(i).p(order);
            data_out(i-n_skipped).cp    = data(i).cp(order);
            data_out(i-n_skipped).ref   = data(i).ref;
            data_out(i-n_skipped).filename = data(i).filename;
            
            % Calculate two-surface coordinate based on normalized-rotated-x coord.
            datx = data_out(i-n_skipped).xnrot;
            [maxx,trailing_index] = max(datx);
            x_twosurf = [         datx(1:trailing_index); ...
                         maxx*2 - datx(trailing_index+1:end)];
            x_twosurf = x_twosurf - maxx;
            data_out(i-n_skipped).x_twosurf = x_twosurf;
            
        end
        
        naca_params = naca_params(run_numbers_kept,:);
        naca_params = [naca_params(:,1:4), ...
                       c*ones(length(run_numbers_kept),1), ...
                       naca_params(:,5),];
        
    catch err
        
        delete(hwait);
        rethrow(err);
        
    end

    delete(hwait);

end

    
    
    
    
    
