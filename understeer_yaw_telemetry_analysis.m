% understeer_yaw_sideslip_telemetry.m
load("alameda05022026skidpad.mat"); % Load any logs with steering angle,
% yaw rate, and lat accel channels; rename variables as necessary

%% --- Vehicle Parameters ---
steering_ratio = 4.0;
wheelbase_m = 1.5742;   

%% --- Data Extraction & Interpolation ---
% 1. Alias log names for cleaner scripting
can_log = rotortemp_fullendurance_05_04_2026_1;
imu_log = RaceBoxTrackSessionon04_05_202605_46;

% 2. Extract base time and steering data (primary timeline)
time_can = can_log.Time;
steer_angle_deg = can_log.STEERINGANGLE;

% 3. Extract IMU/GPS data and its original time vector
time_imu = imu_log.Time; 
speed_raw_ms = imu_log.Speed; 
lat_accel_raw_ms2 = imu_log.GForceY .* 9.806; 
yaw_rate_raw_rads = deg2rad(imu_log.GyroZ); 

% 4. Interpolate IMU data to match CANary resolution (main timeline)
speed_ms = interp1(time_imu, speed_raw_ms, time_can, 'linear', 'extrap');
lat_accel_ms2 = interp1(time_imu, lat_accel_raw_ms2, time_can, 'linear', 'extrap');
yaw_rate_meas_rads = interp1(time_imu, yaw_rate_raw_rads, time_can, 'linear', 'extrap'); 

% Prevent divide-by-zero errors in speed-dependent calculations
speed_ms(speed_ms < 0.1) = 0.1; 

%% --- Understeer & Theoretical Yaw Calculations (Steady-State Approximations) ---
% Linear Yaw Rate w/ K_us (compare to pure kinematic)
[understeer_angle_deg, ideal_yaw_rate_rads, K_us_scalar_deg_g] = ...
    understeer_yaw_telemetry_calc(time_can, steer_angle_deg, steering_ratio, speed_ms, wheelbase_m, lat_accel_ms2);

fprintf('Extracted Understeer Gradient (K_us): %.2f deg/g\n', K_us_scalar_deg_g);

% Convert Road Wheel Angle to Radians for calculation
road_wheel_rad = deg2rad(steer_angle_deg ./ steering_ratio);

% Pure Kinematic Yaw Rate
% Ideal according to only tire angle and wheelbase
yaw_rate_kinematic_rads = (speed_ms .* tan(road_wheel_rad)) ./ wheelbase_m;

% Lateral Acceleration Yaw Rate
% Steady state approximation (r = a_y / v_x)
yaw_rate_lataccel_rads = lat_accel_ms2 ./ speed_ms;

% Sideslip rate of change (transient)
% B = a_y / v_x - measured yaw rate
sideslip_rate_rads = yaw_rate_lataccel_rads - yaw_rate_meas_rads;

%% --- Plots ---

% Subplot 1: Yaw Rate Overlays
subplot(3,1,1);
plot(time_can, yaw_rate_meas_rads, 'LineWidth', 1.5, 'DisplayName', 'Measured');
hold on;
plot(time_can, yaw_rate_kinematic_rads, 'LineWidth', 1.2, 'DisplayName', 'Kinematic (Steering Input)');
plot(time_can, ideal_yaw_rate_rads, 'LineWidth', 1.2, 'DisplayName', 'Theoretical (Linear Tire Model)');
plot(time_can, yaw_rate_lataccel_rads, 'LineWidth', 1.2, 'DisplayName', 'a_y/v_x');

ylabel('Yaw Rate (rad/s)');
title('Yaw Rates');
legend('Location', 'best');
grid on;

% Subplot 2: Instantaneous Understeer Angle
subplot(3,1,2);
plot(time_can, understeer_angle_deg, 'LineWidth', 1.2, 'DisplayName', 'Understeer Angle');
xlabel('Time (s)');
ylabel('Understeer Angle (deg)');
title('Instantaneous Vehicle Understeer Angle (\alpha_{US})');
grid on;

% Subplot 3: Sideslip Rate of Change (Transient)
subplot(3,1,3);
plot(time_can, sideslip_rate_rads, 'LineWidth', 1.2, 'DisplayName', 'Sideslip Rate of Change');
hold on;
xlabel('Time (s)');
ylabel('\beta^{x} (rad/s)');
title('\beta^{.} Sideslip Rate of Change');
legend('Location', 'best');
grid on;