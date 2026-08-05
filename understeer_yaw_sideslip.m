load("alameda05042026.mat"); % Load logs with steering angle, yaw rate, and lat accel channels

%% --- Vehicle Parameters ---
steering_ratio = 4.0;
wheelbase = 1.5742;   % m

%% --- Data Extraction & Interpolation ---
% 1. Extract base time and steering data (primary timeline)
time = rotortemp_fullendurance_05_04_2026_1.Time;
steer_deg = rotortemp_fullendurance_05_04_2026_1.STEERINGANGLE;
% wheelspeed = (rotortemp_fullendurance_05_04_2026_1.WheelSpeedFrontLeft + rotortemp_fullendurance_05_04_2026_1.WheelSpeedFrontRight) ./ 2 ./ 3.6; % convert km/h to m/s

% 2. Extract IMU/GPS data and its original time vector
time_imu = RaceBoxTrackSessionon04_05_202605_46.Time; 
speed_raw = RaceBoxTrackSessionon04_05_202605_46.Speed; % m/s
lat_accel_raw = RaceBoxTrackSessionon04_05_202605_46.GForceY .* 9.806; % m/s^2
measured_yaw_raw = deg2rad(RaceBoxTrackSessionon04_05_202605_46.GyroZ); % rad/s

% 3. Interpolate IMU data to match CANary resolution
speed_interp = interp1(time_imu, speed_raw, time, 'linear', 'extrap');
lat_accel_interp = interp1(time_imu, lat_accel_raw, time, 'linear', 'extrap');
measured_yaw_interp = interp1(time_imu, measured_yaw_raw, time, 'linear', 'extrap'); % Interpolate IMU yaw to main timeline

% Prevent divide-by-zero errors in speed-dependent calculations
speed_interp(speed_interp < 0.1) = 0.1; 

%% --- Understeer & Theoretical Yaw Calculations (Steady-State Approximations) ---
% Linear Yaw Rate w/ K_us (compare to pure kinematic)
[us_angle_deg, r_theory, K_us_scalar_deg_g, K_us_inst_deg_g] = calcUndersteerAndYawRate(...
    steer_deg, steering_ratio, speed_interp, wheelbase, lat_accel_interp);

fprintf('Extracted Vehicle Understeer Gradient (K_us): %.2f deg/g\n', K_us_scalar_deg_g);

% Convert Road Wheel Angle to Radians for calculation
road_rad = deg2rad(steer_deg ./ steering_ratio);

% METHOD A: Pure Kinematic Yaw Rate
% Ideal according to only tire angle and wheelbase
yaw_rate_kinematic_rad = (speed_interp .* tan(road_rad)) ./ wheelbase;

% METHOD B: Lateral Acceleration Yaw Rate
% Steady State assumption (r = Ay / Vx)
yaw_rate_lataccel_rad = lat_accel_interp ./ speed_interp;

% Sideslip rate of change (transient)
% B = a_y/V_x - true yaw rate
beta_dot_rad = yaw_rate_lataccel_rad - measured_yaw_interp;


%% --- Plots ---

% Subplot 1: Yaw Rate Overlays
subplot(3,1,1);
plot(time, measured_yaw_interp, 'LineWidth', 1.5, 'DisplayName', 'IMU Measured');
hold on;
plot(time, yaw_rate_kinematic_rad, 'LineWidth', 1.2, 'DisplayName', 'Kinematic (Driver Intent)');
plot(time, r_theory, 'LineWidth', 1.2, 'DisplayName', 'Theoretical (Linear Tire Model)');
plot(time, yaw_rate_lataccel_rad, 'LineWidth', 1.2, 'DisplayName', 'a_y/V_x Approx)');

ylabel('Yaw Rate (rad/s)');
title('Yaw Rates');
legend('Location', 'best');
grid on;

% Subplot 2: Calculated Understeer Angle
subplot(3,1,2);
plot(time, us_angle_deg, 'LineWidth', 1.2, 'DisplayName', 'Understeer Angle');
xlabel('Time (s)');
ylabel('Understeer Angle (deg)');
title('Instantaneous Vehicle Understeer Angle (\alpha_{US})');
grid on;

% Subplot 3: Beta Dot (Transient)
subplot(3,1,3);
plot(time, beta_dot_rad, 'LineWidth', 1.2, 'Color', '#D95319', 'DisplayName', 'Beta Dot');
hold on;
yline(0, 'LineWidth', 1.5); % The target steady-state line (Zero)
xlabel('Time (s)');
ylabel('Beta Dot (rad/s)');
title('Sideslip Rate of Change (\beta^{.})');
legend('Location', 'best');
grid on;


%% --- Helper Function ---
function [us_angle_deg, r_theory, K_us_scalar_deg_g, K_us_inst_deg_g] = calcUndersteerAndYawRate(delta_sw, sr, vx, L, ay)
    g = 9.80665; % Gravity (m/s^2)

    % 1. Convert Road Wheel Angle to Radians
    delta_wheel_deg = delta_sw ./ sr;
    delta_wheel_rad = deg2rad(delta_wheel_deg);

    % 2. Speed threshold protection
    vx_safe = vx;
    vx_safe(abs(vx_safe) < 0.5) = NaN; 

    % 3. Calculate Ackermann Steer Angle using Lateral Acceleration (Rad)
    delta_ack_rad = (L .* ay) ./ (vx_safe .^ 2);

    % 4. Calculate Instantaneous Understeer Angle
    us_angle_rad = delta_wheel_rad - delta_ack_rad;
    us_angle_deg = rad2deg(us_angle_rad);

    % 5. Calculate Instantaneous Kus (Deg/g)
    ay_safe = ay;
    ay_safe(abs(ay_safe) < 0.2) = NaN; 
    
    K_us_inst_rad_ms2 = us_angle_rad ./ ay_safe;
    K_us_inst_deg_g = rad2deg(K_us_inst_rad_ms2) .* g;

    % 6. Calculate Overall Scalar Kus via Linear Regression
    % Filters for cornering events where tires are in linear range (1.5 m/s^2 to 5.0 m/s^2)
    valid_idx = (abs(ay) > 1.5) & (abs(ay) < 5.0) & (vx > 3.0); % arbitrary; determine better bounds and/or method
    
    if sum(valid_idx) > 20
        P = polyfit(ay(valid_idx), us_angle_rad(valid_idx), 1);
        K_us_scalar_SI = P(1); 
    else
        warning('Insufficient steady-state cornering data for regression. Defaulting Kus to 0 (Neutral Steer).');
        K_us_scalar_SI = 0;
    end
    
    K_us_scalar_deg_g = rad2deg(K_us_scalar_SI) * g;

    %% 7. Calculate Theoretical Yaw Rate (rad/s)
    r_theory = (vx .* delta_wheel_rad) ./ (L + (K_us_scalar_SI .* (vx .^ 2)));
end