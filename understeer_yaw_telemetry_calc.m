function [understeer_angle_deg, ideal_yaw_rate_rads, K_us_scalar_deg_g] = understeer_yaw_telemetry_calc(time_s, steer_wheel_deg, steering_ratio, speed_ms, wheelbase_m, lat_accel_ms2)
% understeer_yaw_telemetry_calc calculates understeer gradients and theoretical yaw rate.
% Requires a time_s vector for exact numerical differentiation.

    g_ms2 = 9.80665;

    % 1. Convert Road Wheel Angle to Radians
    road_wheel_deg = steer_wheel_deg ./ steering_ratio;
    road_wheel_rad = deg2rad(road_wheel_deg);

    % 2. Speed threshold protection
    speed_safe_ms = speed_ms;
    speed_safe_ms(abs(speed_safe_ms) < 5.0) = NaN; % Below 5 m/s, kinematics dominate over dynamics

    % 3. Calculate Ackermann Steer Angle (Rad)
    ackermann_angle_rad = (wheelbase_m .* lat_accel_ms2) ./ (speed_safe_ms .^ 2);

    % 4. Calculate Understeer Angle
    understeer_angle_rad = road_wheel_rad - ackermann_angle_rad;
    understeer_angle_deg = rad2deg(understeer_angle_rad);    

    % 5. Calculate scalar K_us via linear regression
    % Calculate signal rates using the exact time vector to find steady-state conditions
    steer_velocity_deg_s = gradient(steer_wheel_deg, time_s);
    lat_jerk_ms3 = gradient(lat_accel_ms2, time_s);

    % --- Linear tire window filters--EDIT ---
    cond_speed = speed_ms > 6; % filter out kinematic-dominant states
    cond_ay_max = abs(lat_accel_ms2) < 7.7; % 0.5 * skidpad lat. accel from comp (9.125 m radius)
    cond_ss_steer = abs(steer_velocity_deg_s) < 40.0; 
    cond_ss_ay = abs(lat_jerk_ms3) < 2.0; 
    
    ss_pts = cond_speed & cond_ay_max & cond_ss_steer & cond_ss_ay;
    
    if sum(ss_pts) > 50 
        robust_fit_coeffs = robustfit(lat_accel_ms2(ss_pts), understeer_angle_rad(ss_pts));
        K_us_scalar_SI = robust_fit_coeffs(2);
    else
        warning('Insufficient steady-state linear cornering data. Defaulting K_us to B26 target (0.15).');
        K_us_scalar_SI = 0.15;
    end
    
    K_us_scalar_deg_g = rad2deg(K_us_scalar_SI) * g_ms2;

    % 6. Calculate Theoretical Yaw Rate (rad/s)
    ideal_yaw_rate_rads = (speed_ms .* road_wheel_rad) ./ (wheelbase_m + (K_us_scalar_SI .* (speed_ms .^ 2)));
end