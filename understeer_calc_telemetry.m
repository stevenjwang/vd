% understeer_calc_telemetry calculates the scalar understeer gradient K_us
% via linear regression from telemetry logs.

function K_us_scalar_deg_g = understeer_calc_telemetry()

    % Load any skidpad logs with steer angle, v_x, a_x, a_y, and yaw rates;
    % rename variables as necessary
 
    load("alameda05022026skidpad.mat"); 

    % 1. Extract time vectors
    t_racebox = RaceBoxTrackSessionon02_05_202619_52.Time; 
    t_steer   = Alameda_Skidpad0384LLTD_2026_05_02T20_21_47.Time;
    
    % 2. Define common time vector 
    t_common = t_racebox;
    
    % 3. Extract raw variables
    lat_accel_ms2   = RaceBoxTrackSessionon02_05_202619_52.GForceY * 9.80665;
    long_accel_ms2  = RaceBoxTrackSessionon02_05_202619_52.GForceX * 9.80665;
    speed_ms        = RaceBoxTrackSessionon02_05_202619_52.Speed;
    steer_angle_raw = Alameda_Skidpad0384LLTD_2026_05_02T20_21_47.STEERINGANGLE;
    yaw_rate_raw_rads = deg2rad(RaceBoxTrackSessionon02_05_202619_52.GyroZ);
    
    % 4. Interpolate steering angle to match RaceBox timestamps
    steer_angle_deg = interp1(t_steer, steer_angle_raw, t_common, 'linear', 'extrap');
    
    % 5. Define constants
    wheelbase_m = 1.5742;
    steering_ratio = 4;
    g_ms2 = 9.80665;

    % 6. Calculate understeer angle over time (rad)
    road_wheel_rad = deg2rad(steer_angle_deg ./ steering_ratio);
    ideal_angle_rad = (wheelbase_m .* lat_accel_ms2) ./ (speed_ms .^ 2);
    understeer_angle_rad = road_wheel_rad - ideal_angle_rad;

    % 7. Calculate scalar K_us via linear regression
    yaw_accel = gradient(yaw_rate_raw_rads, t_common);

    % --- Linear tire window and steady-state filters ---
    cond_speed = speed_ms > 6; % filter kinematic-dominant states
    cond_ay_min = abs(lat_accel_ms2) > 4; % filter straights
    cond_ay_max = abs(lat_accel_ms2) < 10; % tire linear range
    cond_ax_max = abs(long_accel_ms2) < 2; % filter combined slip
    cond_yaw_accel_max = abs(yaw_accel) < 0.5; % filter turn-in and exit
    
    ss_pts = cond_speed & cond_ax_max & cond_ay_min & cond_ay_max & cond_yaw_accel_max;
    
    if sum(ss_pts) > 50 
        robust_fit_coeffs = robustfit(lat_accel_ms2(ss_pts), understeer_angle_rad(ss_pts));
        K_us_scalar_SI = robust_fit_coeffs(2);
    else
        warning('Insufficient steady-state linear cornering data. Defaulting K_us to B26 target (0.15).');
        K_us_scalar_SI = 0.15;
    end
    
    % Return final scalar value in deg/g
    K_us_scalar_deg_g = rad2deg(K_us_scalar_SI) * g_ms2;

    fprintf('Extracted Understeer Gradient (K_us): %.2f deg/g\n', K_us_scalar_deg_g);
end