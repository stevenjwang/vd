function [x_table_accel, long_vel_guess, long_accel_matrix] = long_accel_sweep(car)
    % Run a sweep of longitudinal accel vs velocity to be used for interpolation
    
    % --- 80/20 OPTIMIZATION ---
    % 50 points provides a smooth enough curve for accurate interpolation
    % while doubling the speed of the sweep compared to 100 points.
    num_points = 50; 
    long_vel_guess = linspace(0.001, car.max_vel, num_points);
    
    % Pre-allocate arrays to prevent memory reallocation lag
    long_accel_matrix = zeros(1, num_points);
    x_cell = cell(num_points, 1); 

    for i = 1:num_points
        if i == 1
            [x_accel, long_accel, long_accel_guess] = max_long_accel(long_vel_guess(i), car);
        else
            % Warm start: Pass the previous optimal state as the new initial guess
            x0 = long_accel_guess;
            [x_accel, long_accel, long_accel_guess] = max_long_accel(long_vel_guess(i), car, x0);
        end
        
        x_cell{i} = x_accel;   
        long_accel_matrix(i) = long_accel;
    end
    
    % Concatenate state matrix instantly
    x_matrix = vertcat(x_cell{:});
    
    % --- BUG FIX ---
    % Truncate the final point (boundary limits at V_max) across ALL arrays
    % so dimensions match perfectly for the interpolant.
    long_vel_guess    = long_vel_guess(1:end-1);
    long_accel_matrix = long_accel_matrix(1:end-1);
    x_matrix          = x_matrix(1:end-1, :); 
    
    [x_table_accel] = generate_table(x_matrix);
    
    % visual check
    % plot(long_vel_guess, long_accel_matrix);
end