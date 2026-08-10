% understeer_yaw_sideslip_telemetry.m
load("alameda05042026.mat"); % Load any time-aligned logs with
% steer angle, v_x, a_x, a_y, and yaw rates; rename logs/channels as necessary

%% --- Vehicle Parameters ---
steering_ratio = 4.0;
wheelbase_m = 1.5742;   

%% --- Data Extraction & Interpolation ---
% 1. Alias log names from mat file
can_log = [rotortemp_fullendurance_05_04_2026_1];
imu_log = [RaceBoxTrackSessionon04_05_202605_46];

% 2. Extract CANary timeline and steering data
time_can = can_log.Time;
steer_angle_deg = can_log.STEERINGANGLE;

% 3. Extract IMU/GPS data and its original time vector
time_imu = imu_log.Time; 
speed_raw_ms = imu_log.Speed; 
lat_accel_raw_ms2 = imu_log.GForceY .* 9.80665; 
long_accel_raw_ms2 = imu_log.GForceX .* 9.80665;
yaw_rate_raw_rads = deg2rad(imu_log.GyroZ); 

% 4. Extract, Filter, and Interpolate IMU data to match CANary resolution 
% Apply moving median to reject spikes
window_size = 3; 

speed_ms = movmedian(speed_raw_ms, window_size);
lat_accel_ms2 = movmedian(lat_accel_raw_ms2, window_size);
long_accel_ms2 = movmedian(long_accel_raw_ms2, window_size);
yaw_rate_rads = movmedian(yaw_rate_raw_rads, window_size);

% Interpolate IMU data
vx_ms = interp1(time_imu, speed_ms, time_can, 'linear', 'extrap');
lat_accel_ms2 = interp1(time_imu, lat_accel_ms2, time_can, 'linear', 'extrap');
long_accel_ms2 = interp1(time_imu, long_accel_ms2, time_can, 'linear', 'extrap');
yaw_rate_meas_deg = rad2deg(interp1(time_imu, yaw_rate_rads, time_can, 'linear', 'extrap'));

% Prevent divide-by-zero errors in speed-dependent calculations
vx_ms(vx_ms < 0.1) = 0.1;

%% --- Understeer & Theoretical Yaw Calculations (Steady-State Approximations) ---

% Convert road wheel angle to radians for calculation
road_wheel_rad = deg2rad(steer_angle_deg ./ steering_ratio);

% Calculate Understeer Gradient (K_us) using standalone function (edit
% params in function)
K_us_scalar_deg_g = understeer_calc_telemetry();

% --- Calculate Understeer Angle ---
% Ideal Ackermann steer angle (rad)
ideal_angle_rad = (wheelbase_m .* lat_accel_ms2) ./ (vx_ms .^ 2);
% Understeer angle
understeer_angle_rad = road_wheel_rad - ideal_angle_rad;
understeer_angle_deg = rad2deg(understeer_angle_rad);

% --- Calculate Theoretical (Ideal) Yaw Rates ---
% Linear Dynamic Yaw Rate (Kinematic with K_us correction)
g_ms2 = 9.80665;
K_us_scalar_SI = deg2rad(K_us_scalar_deg_g) / g_ms2; % Convert back to SI for calculation
ideal_yaw_rate_deg = rad2deg((vx_ms .* road_wheel_rad) ./ (wheelbase_m + (K_us_scalar_SI .* (vx_ms .^ 2))));

% Pure Kinematic Yaw Rate
% Ideal according to only tire angle and wheelbase
yaw_rate_kinematic_deg = rad2deg((vx_ms .* tan(road_wheel_rad)) ./ wheelbase_m);

% Lateral Acceleration Yaw Rate
% Steady state approximation (r = a_y / v_x)
yaw_rate_lataccel_deg = rad2deg(lat_accel_ms2 ./ vx_ms);

% Sideslip rate of change (transient)
% B = a_y / v_x - measured yaw rate
sideslip_rate_deg = rad2deg(yaw_rate_lataccel_deg - yaw_rate_meas_deg);

%% --- Plots ---

% Subplot 1: Yaw Rate Overlays
subplot(4,1,1);
plot(time_can, yaw_rate_meas_deg, 'DisplayName', 'Measured');
hold on;
plot(time_can, yaw_rate_kinematic_deg, 'DisplayName', 'Kinematic (Steering Input)');
plot(time_can, ideal_yaw_rate_deg, 'DisplayName', 'Theoretical (Linear Tire Model)');
plot(time_can, yaw_rate_lataccel_deg, 'DisplayName', 'a_y/v_x');
ylabel('Yaw Rate (deg/s)');
title('Yaw Rates');
legend('Location', 'best');
grid on;


% Subplot 2: Instantaneous Understeer Angle vs Ideal Kinematic Steer Angle
subplot(4,1,2);
plot(time_can, understeer_angle_deg, 'DisplayName', 'Understeer Angle');
xlabel('Time (s)');
ylabel('Understeer Angle (deg)');
title('Instantaneous Understeer Angle (Measured vs Kinematic)');
grid on;

% Subplot 3: Steer Angle
subplot(4,1,3);
plot(time_can, steer_angle_deg, 'DisplayName', 'Steer Angle (deg)');
xlabel('Time (s)');
ylabel('Steer Angle (deg)');
title('Steer Angle');
grid on;


% Subplot 4: Lateral Acceleration
subplot(4,1,4);
plot(time_can, lat_accel_ms2 ./ 9.80665, 'DisplayName', 'a_y');
xlabel('Time (s)');
ylabel('a_y');
title('Lateral Acceleration');
grid on;


% Sideslip Rate of Change (Transient)
figure();
plot(time_can, sideslip_rate_deg, 'DisplayName', 'Sideslip Rate of Change');
hold on;
xlabel('Time (s)');
ylabel('Sideslip Rate of Change (deg/s)');
title('Sideslip Rate of Change');
grid on;
