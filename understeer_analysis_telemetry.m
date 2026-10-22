% Incorporates steady-state K_us from understeer_gradient_telemetry.m, transient criteria 
% (time lags, response times T_90%, overshoot, gain ratios), and sideslip angle.

clear; clc; close all;

load('alameda05042026.mat'); 

steering_ratio = 4.0;    % Steering gear ratio
wheelbase_m    = 1.5742; % Wheelbase L [m]
g_ms2          = 9.80665;% Gravity [m/s^2]
a_dist_m       = 0.80;   % Distance from CG to front axle [m] (Estimated for vehicle geometry)
b_dist_m       = wheelbase_m - a_dist_m; % Distance from CG to rear axle [m]

can_log = rotortemp_fullendurance_05_04_2026_1;
imu_log = RaceBoxTrackSessionon04_05_202605_46;

time_can = can_log.Time;
dt       = mean(diff(time_can));
fs       = 1 / dt; 

steer_angle_deg   = can_log.STEERINGANGLE; % Steering wheel angle \delta_H [deg]
throttle          = can_log.ThrottlePosition;
ws_fl             = can_log.WheelSpeedFrontLeft;
ws_fr             = can_log.WheelSpeedFrontRight;
shock_rl          = can_log.SHOCKRL;
shock_fl          = can_log.SHOCKFL;

time_imu          = imu_log.Time; 
speed_raw_ms      = imu_log.Speed; 
lat_accel_raw_ms2 = imu_log.GForceY .* g_ms2; 
long_accel_raw_ms2= imu_log.GForceX .* g_ms2;
yaw_rate_raw_rads = deg2rad(imu_log.GyroZ); 

% Filter + interpolate
window_size       = 3; 
speed_ms          = movmedian(speed_raw_ms, window_size);
lat_accel_ms2     = movmedian(lat_accel_raw_ms2, window_size);
long_accel_ms2    = movmedian(long_accel_raw_ms2, window_size);
yaw_rate_rads     = movmedian(yaw_rate_raw_rads, window_size);

vx_ms             = interp1(time_imu, speed_ms, time_can, 'linear', 'extrap');
lat_accel_ms2     = interp1(time_imu, lat_accel_ms2, time_can, 'linear', 'extrap');
long_accel_ms2    = interp1(time_imu, long_accel_ms2, time_can, 'linear', 'extrap');
yaw_rate_meas_deg = rad2deg(interp1(time_imu, yaw_rate_rads, time_can, 'linear', 'extrap'));
vx_ms(vx_ms < 0.5) = 0.5; % Prevent zero-division at standstill

% Convert Angles to Radians for Dynamics Formulations
road_wheel_deg    = steer_angle_deg ./ steering_ratio;
road_wheel_rad    = deg2rad(road_wheel_deg);
yaw_rate_meas_rad = deg2rad(yaw_rate_meas_deg);

% Extract steady-state understeer gradient K_us from skidpad
K_us_scalar_deg_g = understeer_gradient_telemetry(); 
K_us_scalar_SI    = deg2rad(K_us_scalar_deg_g) / g_ms2; % [rad / (m/s^2)]

% Yaw Rate Calculations:
% 1. Pure Kinematic (Low-speed geometric ackermann assumption)
yaw_rate_kinematic_deg = rad2deg((vx_ms .* tan(road_wheel_rad)) ./ wheelbase_m);

% 2. Linear Dynamic Theoretical (Kinematic + K_us Understeer Correction)
ideal_yaw_rate_deg = rad2deg((vx_ms .* road_wheel_rad) ./ (wheelbase_m + (K_us_scalar_SI .* (vx_ms .^ 2))));

% 3. Lateral Acceleration Equivalent Yaw Rate (Steady-State Trajectory r = a_y / v_x)
yaw_rate_lataccel_deg = rad2deg(lat_accel_ms2 ./ vx_ms);

% Kinetic Sideslip Rate: \dot{\beta} = (a_y / v_x) - r [rad/s]
sideslip_rate_rads = (lat_accel_ms2 ./ vx_ms) - yaw_rate_meas_rad;
sideslip_rate_degs = rad2deg(sideslip_rate_rads);

% High-Pass Drift-Corrected Numerical Integration for Continuous Sideslip Angle \beta(t)
% Discrete-time high-pass filter (cutoff ~ 0.05 Hz) removes accelerometer/gyro DC bias integration drift
[b_hp, a_hp] = butter(1, 0.05 / (fs / 2), 'high');
raw_sideslip_rad = cumtrapz(time_can, sideslip_rate_rads);

% Zero-reset condition on straightaways (|a_y| < 0.5 m/s^2 & |\delta_H| < 1.0 deg)
straight_idx = (abs(lat_accel_ms2) < 0.5) & (abs(steer_angle_deg) < 1.0);
raw_sideslip_rad(straight_idx) = 0;
sideslip_angle_deg = rad2deg(filtfilt(b_hp, a_hp, raw_sideslip_rad));

% Steady-state baseline sideslip angle \beta_ss for comparison
sideslip_ss_deg = rad2deg((b_dist_m / wheelbase_m) .* road_wheel_rad - (a_dist_m .* lat_accel_ms2) ./ (vx_ms.^2));

% --- A. Sensitivity Gains / Ratios ---
safe_steer_deg = steer_angle_deg;
safe_steer_deg(abs(safe_steer_deg) < 0.5) = NaN; % Deadband filter to prevent divide-by-zero

% Lateral Acceleration divided by Steering Wheel Angle [ (m/s^2) / deg ] and [ g / deg ]
lat_accel_gain_ms2_deg = lat_accel_ms2 ./ safe_steer_deg;
lat_accel_gain_g_deg   = (lat_accel_ms2 ./ g_ms2) ./ safe_steer_deg;

% Yaw Velocity divided by Steering Wheel Angle [ (deg/s) / deg ]
yaw_velocity_gain = yaw_rate_meas_deg ./ safe_steer_deg;

% --- B. Time Lags via Global Cross-Correlation ---
[c_yaw, lags_yaw] = xcorr(yaw_rate_meas_deg - mean(yaw_rate_meas_deg), ...
                          steer_angle_deg - mean(steer_angle_deg), 'normalized');
[~, idx_yaw] = max(c_yaw);
time_lag_steer_yaw_s = lags_yaw(idx_yaw) * dt;

[c_lat, lags_lat] = xcorr(lat_accel_ms2 - mean(lat_accel_ms2), ...
                          steer_angle_deg - mean(steer_angle_deg), 'normalized');
[~, idx_lat] = max(c_lat);
time_lag_steer_lat_s = lags_lat(idx_lat) * dt;

[c_yaw_lat, lags_yaw_lat] = xcorr(lat_accel_ms2 - mean(lat_accel_ms2), ...
                                  yaw_rate_meas_deg - mean(yaw_rate_meas_deg), 'normalized');
[~, idx_yaw_lat] = max(c_yaw_lat);
time_lag_yaw_lat_s = lags_yaw_lat(idx_yaw_lat) * dt;

fprintf('--- Global Time Lags ---\n');
fprintf('Steering Wheel Angle -> Yaw Velocity Time Lag:      %.3f s (%.1f ms)\n', time_lag_steer_yaw_s, time_lag_steer_yaw_s*1000);
fprintf('Steering Wheel Angle -> Lateral Acceleration Time Lag: %.3f s (%.1f ms)\n', time_lag_steer_lat_s, time_lag_steer_lat_s*1000);
fprintf('Yaw Velocity -> Lateral Acceleration Time Lag:      %.3f s (%.1f ms)\n', time_lag_yaw_lat_s, time_lag_yaw_lat_s*1000);

% --- C. Transient Event Extraction (Step Response Criteria: T_90% & Overshoot) ---
% Detect sharp turn-in transient events from continuous telemetry
d_steer_dt = gradient(steer_angle_deg, dt);
steer_transient_peaks = find(abs(d_steer_dt) > 100 & abs(steer_angle_deg) > 15);

% Select the primary isolated transient turn-in maneuver for step response evaluation
if ~isempty(steer_transient_peaks)
    event_idx = steer_transient_peaks(1);
    win_back = round(0.5 / dt); win_fwd = round(2.5 / dt);
    idx_start = max(1, event_idx - win_back);
    idx_end   = min(length(time_can), event_idx + win_fwd);
    
    t_ev = time_can(idx_start:idx_end) - time_can(idx_start);
    st_ev = steer_angle_deg(idx_start:idx_end);
    yr_ev = yaw_rate_meas_deg(idx_start:idx_end);
    ay_ev = lat_accel_ms2(idx_start:idx_end);
    
    % Steering 50% input level timestamp (t_0 per Fig 1)
    st_amp = max(abs(st_ev)) - min(abs(st_ev));
    st_50_val = min(abs(st_ev)) + 0.5 * st_amp;
    idx_t0 = find(abs(st_ev) >= st_50_val, 1, 'first');
    t_0 = t_ev(idx_t0);
    
    % Steady-state levels (averaged over tail of maneuver window)
    tail_pts = round(0.5 / dt);
    yr_ss = mean(yr_ev(end-tail_pts:end));
    ay_ss = mean(ay_ev(end-tail_pts:end));
    
    % Response Times (Time to reach 90% of steady-state level minus t_0)
    idx_yr_90 = find(abs(yr_ev) >= 0.9 * abs(yr_ss), 1, 'first');
    idx_ay_90 = find(abs(ay_ev) >= 0.9 * abs(ay_ss), 1, 'first');
    
    T_90_yaw = t_ev(idx_yr_90) - t_0;
    T_90_lat = t_ev(idx_ay_90) - t_0;
    
    % Overshoot Values: (Peak - SteadyState) / SteadyState * 100%
    yr_peak = max(abs(yr_ev));
    ay_peak = max(abs(ay_ev));
    
    overshoot_yaw_pct = max(0, (yr_peak - abs(yr_ss)) / abs(yr_ss)) * 100;
    overshoot_lat_pct = max(0, (ay_peak - abs(ay_ss)) / abs(ay_ss)) * 100;
    
    fprintf('\n--- Isolated Step Transient Criteria ---\n');
    fprintf('Yaw Velocity Response Time (T_90%%):        %.3f s (%.1f ms)\n', T_90_yaw, T_90_yaw*1000);
    fprintf('Lateral Acceleration Response Time (T_90%%): %.3f s (%.1f ms)\n', T_90_lat, T_90_lat*1000);
    fprintf('Yaw Velocity Overshoot:                    %.1f %%\n', overshoot_yaw_pct);
    fprintf('Lateral Acceleration Overshoot:            %.1f %%\n', overshoot_lat_pct);
else
    t_ev = []; T_90_yaw = NaN; T_90_lat = NaN; overshoot_yaw_pct = NaN; overshoot_lat_pct = NaN;
end

% --- Figure 1: Time-Domain Overlays & Sideslip Dynamics ---
figure('Name', 'Dynamic Responses & Sideslip');
subplot(4,1,1);
plot(time_can, steer_angle_deg);
ylabel('Steer Angle \delta_H (deg)');
title('Steering Wheel Input (\delta_H)');
grid on; xlim([time_can(1) time_can(end)]);

subplot(4,1,2);
plot(time_can, yaw_rate_meas_deg, 'DisplayName', 'Measured'); hold on;
plot(time_can, yaw_rate_kinematic_deg, 'DisplayName', 'Kinematic');
plot(time_can, ideal_yaw_rate_deg, 'DisplayName', 'K_us corrected');
plot(time_can, yaw_rate_lataccel_deg, 'DisplayName', 'Trajectory');
ylabel('Yaw Velocity (deg/s)');
title('Yaw Velocity Overlays');
legend('Location', 'best');
grid on; xlim([time_can(1) time_can(end)]);

subplot(4,1,3);
plot(time_can, lat_accel_ms2 ./ g_ms2);
ylabel('Lat Accel a_y (g)');
title('Lateral Acceleration');
grid on; xlim([time_can(1) time_can(end)]);

subplot(4,1,4);
plot(time_can, sideslip_angle_deg, 'DisplayName', 'Integrated Sideslip Rate (HP Filtered)'); hold on;
plot(time_can, sideslip_ss_deg, 'DisplayName', 'Steady-State Baseline');
plot(time_can, sideslip_rate_degs, 'DisplayName', 'Sideslip Rate (deg/s)');
ylabel('Sideslip Angle'); xlabel('Time (s)');
title('Sideslip Angle vs Sideslip Rate');
legend('Location', 'best', 'Orientation', 'horizontal');
grid on; xlim([time_can(1) time_can(end)]);

% --- Figure 2: Dedicated 3D G-G-V Performance Envelope ---
figure('Name', '3D G-G-V Performance Envelope');
scatter3(long_accel_ms2 ./ g_ms2, lat_accel_ms2 ./ g_ms2, vx_ms .* 3.6, 10, vx_ms .* 3.6, 'filled');
xlabel('Longitudinal Acceleration a_x (g)');
ylabel('Lateral Acceleration a_y (g)');
zlabel('Speed v_x (km/h)');
title('3D G-G-V Diagram (Vehicle Performance Envelope)');
colormap jet; 
cb = colorbar; 
title(cb, 'Speed (km/h)');
grid on;
view(-37.5, 30); 
axis tight;

% --- Figure 3: Sensitivity Gains, Time Lag Spectra, & Transient Step Analysis ---
figure('Name', 'Gains, Time Lags & Transient Criteria');
subplot(2,2,1);
scatter(abs(lat_accel_ms2 ./ g_ms2), abs(yaw_velocity_gain), 8, vx_ms, 'filled');
colormap parula; cb1 = colorbar; title(cb1, 'v_x (m/s)');
xlabel('Lateral Acceleration |a_y| (g)');
ylabel('Yaw Gain r / \delta_H [(deg/s)/deg]');
title('Yaw Velocity Gain vs. Lat Accel');
grid on; ylim([0 2]);

subplot(2,2,2);
scatter(abs(lat_accel_ms2 ./ g_ms2), abs(lat_accel_gain_g_deg), 8, vx_ms, 'filled');
cb2 = colorbar; title(cb2, 'v_x (m/s)');
xlabel('Lateral Acceleration |a_y| (g)');
ylabel('Lat Accel Gain (a_y/g) / \delta_H [g/deg]');
title('Lateral Acceleration Gain vs. Lat Accel');
grid on; ylim([0 0.05]);

subplot(2,2,3);
plot(lags_yaw * dt, c_yaw, 'b-', 'LineWidth', 1.3, 'DisplayName', '\delta_H \rightarrow Yaw Rate'); hold on;
plot(lags_lat * dt, c_lat, 'r--', 'LineWidth', 1.3, 'DisplayName', '\delta_H \rightarrow Lat Accel');
xline(time_lag_steer_yaw_s, 'b:', sprintf('T_{lag,r} = %.1f ms', time_lag_steer_yaw_s*1000));
xline(time_lag_steer_lat_s, 'r:', sprintf('T_{lag,a_y} = %.1f ms', time_lag_steer_lat_s*1000));
xlabel('Time Lag \tau (s)'); ylabel('Normalized Cross-Correlation');
title('Cross-Correlation Time Lags');
legend('Location', 'best'); grid on; xlim([-0.5 0.5]);

subplot(2,2,4);
if ~isempty(steer_transient_peaks)
    plot(t_ev, abs(st_ev) ./ max(abs(st_ev)), 'k--', 'DisplayName', 'Normalized Steer Input'); hold on;
    plot(t_ev, abs(yr_ev) ./ max(abs(yr_ev)), 'b-', 'LineWidth', 1.4, 'DisplayName', 'Normalized Yaw Velocity');
    plot(t_ev, abs(ay_ev) ./ max(abs(ay_ev)), 'r-', 'LineWidth', 1.4, 'DisplayName', 'Normalized Lat Accel');
    xline(t_0, 'k:', 't_{0} (50% Steer)');
    xline(t_0 + T_90_yaw, 'b:', sprintf('Yaw T_{90%%}=%.0f ms', T_90_yaw*1000));
    xline(t_0 + T_90_lat, 'r:', sprintf('Lat T_{90%%}=%.0f ms', T_90_lat*1000));
    xlabel('Window Time (s)'); ylabel('Normalized Response');
    title(sprintf('Transient Step Response (OS_{yaw}=%.1f%%, OS_{lat}=%.1f%%)', overshoot_yaw_pct, overshoot_lat_pct));
    legend('Location', 'best'); grid on;
end

%%

% --- Figure 4: Understeer Gradient Curve (K_us vs. Lateral Acceleration) ---
figure('Name', 'Understeer Gradient K_us vs. Lateral Acceleration Curve');
lat_accel_g = lat_accel_ms2 ./ g_ms2;
valid_kus_idx = (abs(lat_accel_g) > 0.15) & (vx_ms > 4.0) & (abs(steer_angle_deg) > 0.5);

% Instantaneous Understeer Gradient Calculation [deg/g]
% K_us = (Steering Road Wheel Angle / Lat Accel in g) - (Ackermann Steering Angle / Lat Accel in g)
% implement reusable gradient_telemetry

steer_per_g = road_wheel_deg(valid_kus_idx) ./ abs(lat_accel_g(valid_kus_idx));
ackermann_per_g = rad2deg((wheelbase_m .* g_ms2) ./ (vx_ms(valid_kus_idx).^2));
kus_instantaneous = steer_per_g - ackermann_per_g;
ay_filtered_g = abs(lat_accel_g(valid_kus_idx));

% Scatter points colored by vehicle speed
scatter(ay_filtered_g, kus_instantaneous, 12, vx_ms(valid_kus_idx), 'filled');
hold on;

% Compute and plot a moving-median trend line for the understeer characteristic curve
[ay_sorted, sort_idx] = sort(ay_filtered_g);
kus_sorted = kus_instantaneous(sort_idx);
window_pts = max(1, round(length(ay_sorted) / 25));
kus_trend = movmedian(kus_sorted, window_pts);
plot(ay_sorted, kus_trend, 'LineWidth', 2.0, 'DisplayName', 'K_{us} Moving Median');

% Reference line for neutral steering behavior (K_us = 0)
yline(0, 'k--', 'LineWidth', 1.0, 'DisplayName', 'Neutral Steer (K_{us} = 0)');

colormap parula;
cb_kus = colorbar;
title(cb_kus, 'v_x (m/s)');
xlabel('Lateral Acceleration |a_y| (g)');
ylabel('Understeer Gradient K_{us} (deg/g)');
title('Steady-State Understeer Gradient Curve (K_{us} vs. |a_y|)');
legend('Location', 'best');
grid on;
ylim([-5 15]);

%%

% --- Figure 5: Diagnostic Chassis Dynamics & Operational Context ---
figure('Name', 'Chassis Diagnostics & Driving Context');
subplot(3,1,1);
plot(time_can, vx_ms .* 3.6, 'b-', 'LineWidth', 1.2, 'DisplayName', 'GPS Speed'); hold on;
plot(time_can, ws_fl, 'r--', 'DisplayName', 'FL Wheel Speed');
plot(time_can, ws_fr, 'g--', 'DisplayName', 'FR Wheel Speed');
ylabel('Speed (km/h)'); title('Speed Profile'); legend('Location', 'best'); grid on; xlim([time_can(1) time_can(end)]);

subplot(3,1,2);
plot(time_can, throttle, 'k-', 'LineWidth', 1.1);
ylabel('Throttle Position (%)'); title('Throttle Input'); grid on; xlim([time_can(1) time_can(end)]);

subplot(3,1,3);
plot(time_can, shock_rl, 'r-', 'DisplayName', 'Shock RL'); hold on;
plot(time_can, shock_fl, 'b-', 'DisplayName', 'Shock FL');
ylabel('Displacement (mm)'); xlabel('Time (s)'); title('Damper Travel'); legend('Location', 'best'); grid on; xlim([time_can(1) time_can(end)]);

% --- Figure 6: K_us vs. Lat Accel & Speed ---
figure('Name', '3D Operational Understeer Map: K_us vs. Lat Accel & Speed');
scatter3(abs(lat_accel_g(valid_kus_idx)), vx_ms(valid_kus_idx) .* 3.6, kus_instantaneous, 14, kus_instantaneous, 'filled');
xlabel('Lateral Acceleration |a_y| (g)');
ylabel('Vehicle Speed v_x (km/h)');
zlabel('Understeer Gradient K_{us} (deg/g)');
title('3D Operational Understeer Map (K_{us} vs. |a_y| & v_x)');
colormap jet; 
cb_3d = colorbar; 
title(cb_3d, 'K_{us} (deg/g)');
grid on; 
view(-40, 25); 
axis tight;

figure('Name', 'Sideslip vs. Lateral Acceleration');

% Subplot 1: Integrated & Theoretical Sideslip Angle vs. Lat Accel (g)
subplot(3,1,1);
valid_slip_idx = (abs(lat_accel_g) > 0.05) & (vx_ms > 2.0);
scatter(abs(lat_accel_g(valid_slip_idx)), abs(sideslip_angle_deg(valid_slip_idx)), 10, vx_ms(valid_slip_idx), 'filled');
hold on;
scatter(abs(lat_accel_g(valid_slip_idx)), abs(sideslip_ss_deg(valid_slip_idx)), 6, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.4);

% Median Trend line for Integrated Sideslip Angle
[ay_slip_sorted, slip_sort_idx] = sort(abs(lat_accel_g(valid_slip_idx)));
beta_sorted = abs(sideslip_angle_deg(valid_slip_idx));
beta_sorted = beta_sorted(slip_sort_idx);
window_slip_pts = max(1, round(length(ay_slip_sorted) / 30));
beta_trend = movmedian(beta_sorted, window_slip_pts);

plot(ay_slip_sorted, beta_trend, 'r-', 'LineWidth', 2.0, 'DisplayName', '|\beta| Dynamic Trend');
colormap parula;
cb_s1 = colorbar; title(cb_s1, 'v_x (m/s)');
xlabel('Lateral Acceleration |a_y| (g)');
ylabel('Sideslip Angle |\beta| (deg)');
title('Sideslip Angle vs. Lateral Acceleration');
legend('|\beta| Dynamic', '|\beta_{ss}| Baseline', '|\beta| Trend', 'Location', 'northwest');
grid on;

% Subplot 2: Kinetic Sideslip Rate vs. Lat Accel (g)
subplot(3,1,2);
scatter(abs(lat_accel_g(valid_slip_idx)), abs(sideslip_rate_degs(valid_slip_idx)), 10, vx_ms(valid_slip_idx), 'filled');
hold on;
beta_rate_sorted = abs(sideslip_rate_degs(valid_slip_idx));
beta_rate_sorted = beta_rate_sorted(slip_sort_idx);
beta_rate_trend = movmedian(beta_rate_sorted, window_slip_pts);

plot(ay_slip_sorted, beta_rate_trend, 'm-', 'LineWidth', 2.0, 'DisplayName', '|\dot{\beta}| Trend');
cb_s2 = colorbar; title(cb_s2, 'v_x (m/s)');
xlabel('Lateral Acceleration |a_y| (g)');
ylabel('Sideslip Rate |\dot{\beta}| (deg/s)');
title('Sideslip Rate vs. Lateral Acceleration');
legend('Location', 'northwest');
grid on;

% Subplot 3: Sideslip Sensitivity Gain (\beta / a_y) vs. Lat Accel (g)
subplot(3,1,3);
beta_gain_deg_g = abs(sideslip_angle_deg(valid_slip_idx)) ./ abs(lat_accel_g(valid_slip_idx));
% Filter extreme spikes from zero-crossing divisions
valid_gain_mask = beta_gain_deg_g < 20;

scatter(abs(lat_accel_g(valid_slip_idx(valid_gain_mask))), beta_gain_deg_g(valid_gain_mask), 10, vx_ms(valid_slip_idx(valid_gain_mask)), 'filled');
hold on;

[ay_gain_sorted, gain_sort_idx] = sort(abs(lat_accel_g(valid_slip_idx(valid_gain_mask))));
gain_sorted = beta_gain_deg_g(valid_gain_mask);
gain_sorted = gain_sorted(gain_sort_idx);
gain_trend = movmedian(gain_sorted, max(1, round(length(ay_gain_sorted) / 30)));

plot(ay_gain_sorted, gain_trend, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Gain Trend');
cb_s3 = colorbar; title(cb_s3, 'v_x (m/s)');
xlabel('Lateral Acceleration |a_y| (g)');
ylabel('Sideslip Gain |\beta / a_y| (deg/g)');
title('Sideslip Sensitivity Gain vs. Lateral Acceleration');
legend('Location', 'northeast');
grid on; ylim([0 15]);