function [understeer_angle_deg, ideal_yaw_rate_rads, K_us_scalar_deg_g, K_us_inst_deg_g] = understeer_yaw_telemetry_calc(steer_wheel_deg, steering_ratio, speed_ms, wheelbase_m, lat_accel_ms2)
% understeer_yaw_calc calculates instantaneous and scalar understeer gradients and theoretical yaw rate

    g_ms2 = 9.80665;

    % 1. Convert Road Wheel Angle to Radians
    road_wheel_deg = steer_wheel_deg ./ steering_ratio;
    road_wheel_rad = deg2rad(road_wheel_deg);

    % 2. Speed threshold protection
    speed_safe_ms = speed_ms;
    speed_safe_ms(abs(speed_safe_ms) < 0.5) = NaN; 

    % 3. Calculate Ackermann Steer Angle using Lateral Acceleration (Rad)
    ackermann_angle_rad = (wheelbase_m .* lat_accel_ms2) ./ (speed_safe_ms .^ 2);

    % 4. Calculate Instantaneous Understeer Angle
    understeer_angle_rad = road_wheel_rad - ackermann_angle_rad;
    understeer_angle_deg = rad2deg(understeer_angle_rad);

    % % 5. Calculate Instantaneous K_us (Deg/g)
    % % INSTANTANEOUS: d(understeer angle) / d(a_y)?
    % lat_accel_safe_ms2 = lat_accel_ms2;
    % lat_accel_safe_ms2(abs(lat_accel_safe_ms2) < 0.2) = NaN; 
    % 
    % K_us_inst_rad_ms2 = understeer_angle_rad ./ lat_accel_safe_ms2;
    % K_us_inst_deg_g = rad2deg(K_us_inst_rad_ms2) .* g_ms2;

    

    % 5. Calculate Instantaneous K_us (Deg/g)
    % Calculated as d(understeer angle) / d(a_y) using discrete gradients
    
    % Get the rate of change between data points
    d_understeer_rad = gradient(understeer_angle_rad);
    d_lat_accel = gradient(lat_accel_ms2);
    
    % Prevent division by near-zero changes in lateral accel (avoids noise/infinity spikes)
    d_lat_accel_safe = d_lat_accel;
    d_lat_accel_safe(abs(d_lat_accel_safe) < 1e-3) = NaN; 
    
    % Calculate the tangent slope
    K_us_inst_rad_ms2 = d_understeer_rad ./ d_lat_accel_safe;
    
    % Convert rad/(m/s^2) to deg/g
    K_us_inst_deg_g = rad2deg(K_us_inst_rad_ms2) .* g_ms2;



    % 6. Calculate Overall Scalar K_us via Linear Regression
    % CHANGE filters for cornering events where tires are in linear range (1.5 m/s^2 to 5.0 m/s^2)
    % something like abs(yaw_rate) < 2
    valid_idx = (abs(lat_accel_ms2) > 1.5) & (abs(lat_accel_ms2) < 5.0) & (speed_ms > 3.0); 
    
    if sum(valid_idx) > 20
        poly_fit_coeffs = polyfit(lat_accel_ms2(valid_idx), understeer_angle_rad(valid_idx), 1);
        K_us_scalar_SI = poly_fit_coeffs(1); 
    else
        warning('Insufficient steady-state cornering data for regression. Defaulting K_us to 0 (Neutral Steer).');
        K_us_scalar_SI = 0;
    end
    
    K_us_scalar_deg_g = rad2deg(K_us_scalar_SI) * g_ms2;

    % 7. Calculate Theoretical Yaw Rate (rad/s)
    ideal_yaw_rate_rads = (speed_ms .* road_wheel_rad) ./ (wheelbase_m + (K_us_scalar_SI .* (speed_ms .^ 2)));
end