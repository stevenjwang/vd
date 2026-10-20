function K_us_scalar_deg_g = understeer_gradient_telemetry()
% Calculates scalar understeer gradient K_us from any time-aligned log containing
% steer angle, v_x, a_x, a_y, and yaw rate; rename variables as necessary

load("alameda05022026skidpad.mat"); % Skidpad gives best K_us conditions;
% update as needed

% Alias log names for cleaner scripting
can_log = [Alameda_Skidpad0384LLTD_2026_05_02T20_21_47];
imu_log = [RaceBoxTrackSessionon02_05_202619_52];

% Ensure logs are time-aligned
% can_log = readtable(Alameda_Skidpad0384LLTD_2026_05_02T20_21_47);
% imu_log = readtable(RaceBoxTrackSessionon02_05_202619_52);

% 1. Extract time vectors
t_racebox = imu_log.Time;
t_canary = can_log.Time;

% 2. Filtered variables
window_size = 3;
lat_accel_ms2     = movmedian(imu_log.GForceY * 9.80665, window_size); % May need to update channel names
long_accel_ms2    = movmedian(imu_log.GForceX * 9.80665, window_size);
vx_ms             = movmedian(imu_log.Speed, window_size);
yaw_rate_rads = movmedian(deg2rad(imu_log.GyroZ), window_size);
steer_angle       = can_log.STEERINGANGLE;

% 3. Interpolate IMU channels to CANary
lat_accel_ms2     = interp1(t_racebox, lat_accel_ms2, t_canary, 'linear', 'extrap');
long_accel_ms2    = interp1(t_racebox, long_accel_ms2, t_canary, 'linear', 'extrap');
vx_ms             = interp1(t_racebox, vx_ms, t_canary, 'linear', 'extrap');
yaw_rate_rads_interp = interp1(t_racebox, yaw_rate_rads, t_canary, 'linear', 'extrap');

% 4. Constants
wheelbase_m = 1.5742;
steering_ratio = 4;
g_ms2 = 9.80665;

% 5. Understeer angle (rad)
road_wheel_rad = deg2rad(steer_angle ./ steering_ratio);
ideal_angle_rad = (wheelbase_m .* lat_accel_ms2) ./ (vx_ms .^ 2);
understeer_angle_rad = road_wheel_rad - ideal_angle_rad;

% 6. Filters: pretty magical numbers to approximate steady-state
yaw_accel   = sgolayfilt(gradient(yaw_rate_rads_interp, t_canary), 3, 11);
steer_rate  = sgolayfilt(gradient(deg2rad(steer_angle), t_canary), 3, 11);

cond_speed = vx_ms > 6; % filter out kinematic-dominant states
cond_ay_min = abs(lat_accel_ms2) > 0.6; % 4 = lower bound of a_y grouping near a_x = 0 from endurance g-g plot
cond_ay_max = abs(lat_accel_ms2) < 10; % conservative LC0 linear range
cond_ax_max = abs(long_accel_ms2) < 0.3; % 2 = a_x ~= 0 grouping from endurance g-g plot
cond_yaw_accel  = abs(yaw_accel) < 0.15;
cond_steer_rate = abs(steer_rate) < 0.1;

ss_pts = cond_speed & cond_ax_max & cond_ay_min & cond_ay_max & cond_yaw_accel & cond_steer_rate;

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