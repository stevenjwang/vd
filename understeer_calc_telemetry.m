function K_us_scalar_deg_g = understeer_calc_telemetry()
% Calculates scalar understeer gradient K_us from any log containing steer angle,
% v_x, a_x, a_y, and yaw rate; rename variables as necessary

load("alameda05022026skidpad.mat"); % Skidpad gives best K_us conditions

% Alias log names for cleaner scripting
can_log = [Alameda_Skidpad0384LLTD_2026_05_02T20_21_47];
imu_log = [RaceBoxTrackSessionon02_05_202619_52];

% 1. Extract time vectors
t_racebox = imu_log.Time;
t_canary = can_log.Time;

% 2. Common time vector
t_common = t_racebox;

% 3. Raw variables (assumes same field names as original logs)
lat_accel_ms2   = imu_log.GForceY * 9.80665;
long_accel_ms2  = imu_log.GForceX * 9.80665;
vx_ms        = imu_log.Speed;
steer_angle_raw = can_log.STEERINGANGLE;
yaw_rate_raw_rads = deg2rad(imu_log.GyroZ);

% 4. Interpolate steering angle
steer_angle_deg = interp1(t_canary, steer_angle_raw, t_common, 'linear', 'extrap');

% 5. Constants
wheelbase_m = 1.5742;
steering_ratio = 4;
g_ms2 = 9.80665;

% 6. Understeer angle (rad)
road_wheel_rad = deg2rad(steer_angle_deg ./ steering_ratio);
ideal_angle_rad = (wheelbase_m .* lat_accel_ms2) ./ (vx_ms .^ 2);
understeer_angle_rad = road_wheel_rad - ideal_angle_rad;

% 7. K_us via robust linear fit with steady-state filters
yaw_accel = gradient(yaw_rate_raw_rads, t_common);

cond_speed = vx_ms > 6;
cond_ay_min = abs(lat_accel_ms2) > 4;
cond_ay_max = abs(lat_accel_ms2) < 10;
cond_ax_max = abs(long_accel_ms2) < 2;
cond_yaw_accel_max = abs(yaw_accel) < 0.5;

ss_pts = cond_speed & cond_ax_max & cond_ay_min & cond_ay_max & cond_yaw_accel_max;

if sum(ss_pts) > 50
    coeffs = robustfit(lat_accel_ms2(ss_pts), understeer_angle_rad(ss_pts));
    K_us_scalar_SI = coeffs(2);
else
    warning('Insufficient steady-state linear cornering data. Defaulting K_us to B26 target (0.15).');
    K_us_scalar_SI = 0.15;
end

K_us_scalar_deg_g = rad2deg(K_us_scalar_SI) * g_ms2;
fprintf('Extracted Understeer Gradient (K_us): %.2f deg/g\n', K_us_scalar_deg_g);
end
