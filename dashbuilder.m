% understeer_yaw_complement_plots.m

load('alameda05042026.mat');

canlog = rotortemp_fullendurance_05_04_2026_1;
imulog = RaceBoxTrackSessionon04_05_202605_46;
imuspeed = imulog.Speed .* 3.6; % m/s to km/h

figure();

subplot(4,1,1);
plot(imulog.Time,imuspeed,'DisplayName','IMUspeed', XDataSource = 'imulog.Time',YDataSource = 'imuspeed');
linkdata on;
xlabel("Time");
ylabel("Speed (km/h)");
title("Speed vs Time");
hold on;
plot(canlog.Time,canlog.WheelSpeedFrontLeft,'DisplayName','WheelSpeedFrontLeft',XDataSource = 'canlog.Time',YDataSource = 'canlog.WheelSpeedFrontLeft');
hold on;
plot(canlog.Time,canlog.WheelSpeedFrontRight,'DisplayName','WheelSpeedFrontRight',XDataSource = 'canlog.Time',YDataSource = 'canlog.WheelSpeedFrontRight');
linkdata on;
legend("show");
grid on;
hold on;

subplot(4,1,2);
plot(canlog.Time,canlog.ThrottlePosition,XDataSource = 'canlog.Time',YDataSource = 'canlog.ThrottlePosition');
linkdata on;
xlabel("Time");
ylabel("ThrottlePosition");
title("ThrottlePosition vs Time");
grid on;

subplot(4,1,3);
plot(imulog.Time,imulog.GForceX,XDataSource = 'imulog.Time',YDataSource = 'imulog.GForceX');
linkdata on;
xlabel("Time");
ylabel("GForceX");
title("GForceX");
grid on;


subplot(4,1,4);
plot(canlog.Time,canlog.SHOCKRL,'DisplayName','canlog.SHOCKRL');
hold on;
plot(canlog.Time,canlog.SHOCKFL,'DisplayName','canlog.SHOCKFL');
hold off;
linkdata on;
ylabel("Shockpot Length (mm)");
title("Shockpot Lengths");
legend("show");
grid on;