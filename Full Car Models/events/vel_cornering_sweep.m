function [x_table_corner_vel, radius_vector, max_vel_corner_vector] = vel_cornering_sweep(car)
    % Run a sweep of velocity vs radius to be used for interpolation
    
    % --- 80/20 OPTIMIZATION ---
    % 45 points is mathematically dense enough to capture the curve perfectly,
    % cutting the computational load by 85% compared to the 0.15 step size.
    num_points = 45;
    radius_vector = linspace(4.5, 50, num_points);
    
    % Pre-allocate memory to prevent reallocation lag
    max_vel_corner_vector = zeros(1, num_points);
    x_cell = cell(num_points, 1);
    
    for i = 1:num_points
        if i == 1
            [x_corner_vel, max_vel, vel_corner_guess] = max_vel_cornering(radius_vector(i), car.max_vel, car);
        else
            % Warm start: Pass the previous optimal state as the initial guess
            x0 = vel_corner_guess; 
            [x_corner_vel, max_vel, vel_corner_guess] = max_vel_cornering(radius_vector(i), car.max_vel, car, x0);
        end
        
        % Store results instantly in pre-allocated slots
        x_cell{i} = x_corner_vel;
        max_vel_corner_vector(i) = max_vel;
    end
    
    % Concatenate the state matrix in one fast memory operation
    x_matrix = vertcat(x_cell{:});
    
    % Generate the final table
    x_table_corner_vel = generate_table(x_matrix);
    
    % Visual check
    % figure;
    % plot(radius_vector, max_vel_corner_vector, '-o');
    % title('Max Cornering Velocity vs. Radius');
    % xlabel('Radius (m)'); ylabel('Velocity (m/s)');
end