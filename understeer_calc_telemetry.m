function K_us_scalar_deg_g = understeer_calc_telemetry()
% Calculates scalar understeer gradient K_us from any time-aligned log containing
% steer angle, v_x, a_x, a_y, and yaw rate assuming constant radius

% 1. Load data
can_log = readtable('2026-10-04_alameda_skidpad_run1_MR2_long(4)_aero_on.csv');
imu_log = readtable('2026-10-04_alameda_skidpad_run1_MR2_long(4)_aero_on.csv');

% 2. Extract time vectors
t_racebox = imu_log.Time;
t_canary = can_log.Time;

%%

% 3. Filtered variables
window_size = 5;
lat_accel_ms2     = movmedian(imu_log.ACCELY * 9.80665, window_size);
long_accel_ms2    = movmedian(imu_log.ACCELX * 9.80665, window_size);
vx_ms             = movmedian(imu_log.WheelSpeedRearRight, window_size);
yaw_rate_rads     = movmedian(deg2rad(imu_log.ANGRATEY), window_size);
steer_angle       = can_log.STEERINGANGLE;

% 4. Interpolate IMU channels to CANary
lat_accel_ms2        = interp1(t_racebox, lat_accel_ms2, t_canary, 'linear', 'extrap');
long_accel_ms2       = interp1(t_racebox, long_accel_ms2, t_canary, 'linear', 'extrap');
vx_ms                = interp1(t_racebox, vx_ms, t_canary, 'linear', 'extrap');
yaw_rate_rads_interp = interp1(t_racebox, yaw_rate_rads, t_canary, 'linear', 'extrap');
steer_angle_deg      = steer_angle;

% 5. Constants
wheelbase_m = 1.5742;
steering_ratio = 4;
g_ms2 = 9.80665;

% 6. Understeer angle (rad)
road_wheel_rad = deg2rad(steer_angle_deg ./ steering_ratio);
ideal_angle_rad = (wheelbase_m .* lat_accel_ms2) ./ (vx_ms .^ 2);
understeer_angle_rad = road_wheel_rad - ideal_angle_rad;

% 7. Enhanced Steady-State Filtering & Robust Fit
% Derivatives mapped to t_canary
yaw_accel      = gradient(yaw_rate_rads_interp, t_canary);
steering_rate  = abs(gradient(steer_angle_deg, t_canary));
lat_accel_rate = abs(gradient(lat_accel_ms2, t_canary));

% Threshold parameters
min_speed          = 6.0;   % m/s (~21.6 km/h) - filter out low-speed kinematics
max_long_accel     = 2;   % m/s^2 - tight bound for near-constant speed (prevents roll-on/braking transients)
min_lat_accel      = 2.5;   % m/s^2 - lower bound to isolate active loaded cornering
max_lat_accel      = 16.0;  % m/s^2 - raised ceiling to handle high-g aero/tire limits safely
max_yaw_accel      = 0.5;   % rad/s^2 - filters transient yaw rotation
max_steering_rate  = 2.5;   % deg/s - ensures driver hands are steady in the circle
max_lat_jerk       = 1.0;   % m/s^3 - ensures a_y has physically settled into radius

cond_speed      = vx_ms > min_speed;
cond_ax         = abs(long_accel_ms2) < max_long_accel;
cond_ay         = (abs(lat_accel_ms2) > min_lat_accel) & (abs(lat_accel_ms2) < max_lat_accel);
cond_yaw_accel  = abs(yaw_accel) < max_yaw_accel;
cond_steer_rate = steering_rate < max_steering_rate;
cond_lat_jerk   = lat_accel_rate < max_lat_jerk;

% Combined mask
ss_pts = cond_speed & cond_ax & cond_ay & cond_yaw_accel; % & cond_lat_jerk; % & cond_steer_rate;

%%

if sum(ss_pts) > 50
    coeffs = robustfit(lat_accel_ms2(ss_pts), understeer_angle_rad(ss_pts));
    K_us_scalar_SI = coeffs(2);
    K_us_scalar_deg_g = rad2deg(K_us_scalar_SI) * g_ms2;
else
    warning('Insufficient steady-state linear cornering data. Defaulting K_us to B26 target (0.15).');
    K_us_scalar_deg_g = 0.15;
end

fprintf('Extracted Understeer Gradient (K_us): %.2f deg/g\n', K_us_scalar_deg_g);
end