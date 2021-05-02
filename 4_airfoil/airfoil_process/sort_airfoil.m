function [ cp_mat_out, xx] = ...
    sort_airfoil( naca_params, data, x_coords, y_coords, scalar_data)
%%%
% Sorts data points within each entry of the NACA data structure, such that they progress
% around the circumference of the deformed wing. This is done by computing the camber line
% of the deformed airfoil, and identifying nodes above and below it. These two sets of
% nodes are then sorted by increasing and decreasing x-coordinate, respectively.
%%%
%%% Inputs: 
% NACA simulation parameters - n_sim x n_params
% data - pressure on airfoil surface - n_sim x n_points
% x_coords 
% y_coords
% scalar_data - probe pressure for cp calculation. 

% Iterpolate to uniform 200 points:
n_points = 200; 
n_sims = size(data,1); 

xx = linspace(-1, 1, n_points); 
cp_mat_out = zeros(n_points,n_sims);

% to calculate cp: 
rho = 1; 
u_ref = 120; 


for i_sim = 1:n_sims

    m = naca_params(i_sim,1);
    p = naca_params(i_sim,2);
    t = naca_params(i_sim,3);
    a = naca_params(i_sim,4);
    c = naca_params(i_sim,5);

    %%% calculate cp 

    p_ref = scalar_data(i_sim); 
    p_data = data(i_sim,:); 

    cp = (p_data - p_ref) ./ (0.5 * rho * (u_ref)^2);

    %  Rotate surface coordinate pairs back to 0 angle-of-attack.
    arad = a * pi / 180;
    rotmat = [cos(arad), -sin(arad); sin(arad), cos(arad)];
    tmp = rotmat * [x_coords(i_sim,:);y_coords(i_sim,:)];
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


    c_ordered = cp(order); 
    x_1 = xnrot(order); 
    [maxx,trailing_index] = max(x_1);
    x_twosurf = [         x_1(1:trailing_index); ...
                 maxx*2 - x_1(trailing_index+1:end)];
    x_twosurf = x_twosurf - maxx;

    % Some coordinates are repeated. Delete repetitions. 
    [x_twosurf_unique, IA, ~] = unique(x_twosurf);
    c_ordered_unique = c_ordered(IA); 

    cp_mat_out(:,i_sim) = interp1(x_twosurf_unique, c_ordered_unique, xx, 'linear','extrap'); 
end

end

    
    
    
    
    
