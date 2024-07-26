clear
close all

launch_alt = 0;

% from sd card
data = readmatrix("telemetrylog.txt");

lowtime = 18295;
hightime = 152412; %204730;

data(:,45) = (data(:,45) - data(1,45))/10^6;

trim = data(lowtime:hightime,:);
% trim = trim(trim(:,45) < hightime,:);
time = trim(:,47) ./ 1e6;

timeAll = time;
timeStart = find(time > 1032.37);
timeEnd = find(time < 1220.09);
trim = trim(timeStart(1) : timeEnd(end), :);
time = timeAll(timeStart(1) : timeEnd(end), :);
time = time - time(1);

gps_long = trim(:,1);
gps_lat = trim(:,2);
gps_alt = trim(:,3);

gps_vn = trim(:,4);
gps_ve = trim(:,5);
gps_vd = trim(:,6);

gps_sat = trim(:,7);
gps_fix = trim(:,8);

ax = trim(:,9);
ay = trim(:,10);
az = trim(:,11);

h_ax = trim(:,12);
h_ay = trim(:,13);
h_az = trim(:,14);

gx = trim(:,15);
gy = trim(:,16);
gz = trim(:,17);

mx = trim(:,18);
my = trim(:,19);
mz = trim(:,20);

imu_temp = trim(:,21);

baro_alt = trim(:,22);
baro_temp = trim(:,23);
baro_press = trim(:,24);

batt_voltage = trim(:,25);
batt_percent = trim(:,26);

roll = trim(:,27);
pitch = trim(:,28);
yaw = trim(:,29);

q0 = trim(:,30);
q1 = trim(:,31);
q2 = trim(:,32);
q3 = trim(:,33);

pn = trim(:,34);
pe = trim(:,35);
pd = trim(:,36);

vn = trim(:,37);
ve = trim(:,38);
vd = trim(:,39);

an = trim(:,40);
ae = trim(:,41);
ad = trim(:,42);

rssi = trim(:,43);
snr = trim(:,44);

%%
apogee = max(-pd);
apogeeIdx = find(-pd == apogee);
%%

f1 = figure()
hold on 
grid on
grid minor
box on
% plot(time,pn,'LineWidth',2)
% plot(time,pe,'LineWidth',2)
% plot(time,pd,'LineWidth',2)
plot3(pn(1:apogeeIdx), pe(1:apogeeIdx), -pd(1:apogeeIdx),'-','LineWidth',5, 'Color', '#0072BD')
plot3(pn(apogeeIdx + 1 : end), pe(apogeeIdx + 1 : end), -pd(apogeeIdx + 1 : end),'-','LineWidth',5, 'Color', '#D95319')
% xlim([-20,400]);
% ylim([-20,400]);
axis square
% plot(time,-baroalt - 86.79)
legend('Ascent', 'Descent')
hold off

xlabel("North (m)")
ylabel("East (m)")
zlabel("Altitude (m)")

nicePlot(f1)

%%

% Initialize video
myVideo = VideoWriter('myVideoFile', 'MPEG-4'); %open video file
myVideo.FrameRate = 50;  %can adjust this, 5 - 10 works well for me
open(myVideo)

% Set up the left axis
f4 = figure();
box on; grid on;

h1 = animatedline('LineWidth', 2, 'Color', "#0072BD");
% ylabel('Temperature ratio - T_inf/T_e');
% ylim([-0.2 1.0]) % set limits before animation

grid on;
% h = animatedline('MaximumNumPoints', 10);

% Force a 3D view
view(3);
numpoints = 4;

% for k = 1:10:numpoints-99
axis equal    
xlim([-1000 100]); ylim([0 650]); zlim([0 1000]);
xlabel("North (m)", 'Interpreter','latex', 'FontSize',16)
ylabel("East (m)", 'Interpreter','latex', 'FontSize',16)
zlabel("Altitude (m)", 'Interpreter','latex', 'FontSize',16)
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
for k = 1:numpoints:length(time)-numpoints - 1
    xvec = pn(k:k+numpoints-1);
    yvec = pe(k:k+numpoints-1);
    zvec = -pd(k:k+numpoints-1);
    % addpoints(h,xvec,yvec)
    addpoints(h1,xvec, yvec, zvec);
    drawnow
    % pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

%%

f2 = figure();
hold on 
grid on
grid minor
box on
p1 = plot(time, baro_press); % after yyaxis left; set to axis color.
xlabel('Time (s)'); ylabel('Barometric Pressure (Pa)', 'Interpreter','latex');

hold off

% set up the right axis
f3 = figure();
hold on 
grid on
grid minor
box on
p2 = plot(time, baro_temp); % after yyaxis right; set to axis color.
xlabel('Time (s)'); ylabel('Barometric Temperature (K)');
hold off


nicePlot(f2)
nicePlot(f3)
% set(axr, 'FontSize', 16)
% set(axl, 'FontSize', 16)

%%
figure()
hold on 
grid on
grid minor
box on
plot(time,-pd,'-*','LineWidth',2)
hold off

xlabel("Time(s)")
ylabel("Up (m)")

figure()
hold on 
grid on
grid minor
box on
plot(time,an,'-*','LineWidth',2)
plot(time,ae,'-*','LineWidth',2)
plot(time,ad,'-*','LineWidth',2)
hold off
legend('an','ae','ad')

xlabel("Time(s)")
ylabel("Acceleration (m/s)")

figure()
hold on 
grid on
grid minor
box on
plot(time,h_ax,'-*','LineWidth',2)
plot(time,h_ay,'-*','LineWidth',2)
plot(time,h_az,'-*','LineWidth',2)
hold off
legend('ax','ay','az')

xlabel("Time(s)")
ylabel("Acceleration (g)")

figure()
hold on 
grid on
grid minor
box on
plot(time,ax,'-*','LineWidth',2)
plot(time,ay,'-*','LineWidth',2)
plot(time,az,'-*','LineWidth',2)
hold off
legend('ax','ay','az')

xlabel("Time(s)")
ylabel("Acceleration (g)")

figure()
hold on 
grid on
grid minor
box on
p3 = plot(time,yaw,'-*','LineWidth',2);
p1 = plot(time,roll,'-*','LineWidth',2);
p2 = plot(time,pitch,'-*','LineWidth',2);
hold off
legend([p1, p2, p3], 'roll','pitch','yaw')

xlabel("Time(s)")
ylabel("Euler angles")

figure()
hold on 
grid on
grid minor
box on
plot(time,vn,'-*','LineWidth',2)
plot(time,ve,'-*','LineWidth',2)
plot(time,vd,'-*','LineWidth',2)
hold off
legend('vn','ve','vd')

xlabel("Time(s)")
ylabel("Velocity (m/s)")

figure()
hold on 
grid on
grid minor
box on
plot(time,rssi,'-*','LineWidth',2)
hold off

xlabel("Time(s)")
ylabel("RSSI (dB)")

figure()
hold on 
grid on
grid minor
box on
plot(time,snr,'-*','LineWidth',2)
hold off

xlabel("Time(s)")
ylabel("SNR (dB)")

figure()
hold on 
grid on
grid minor
box on
plot(time,(gps_alt-launch_alt)/10^3,'-*','LineWidth',2)
hold off

xlabel("Time(s)")
ylabel("GPS Altitude (m)")

speed = sqrt(vn.^2 + vd.^2 + ve.^2);

figure()
hold on 
grid on
grid minor
box on
plot(time,speed,'-*','LineWidth',2)
hold off
xlabel("Time(s)")
ylabel("Speed (m/s)")

figure()
hold on 
grid on
grid minor
box on
plot(time,-vd,'-*','LineWidth',2)
hold off
xlabel("Time(s)")
ylabel("Descent Velocity (m/s)")

figure()
hold on 
grid on
grid minor
box on
plot(time,sqrt(vn.^2+ve.^2),'-*','LineWidth',2)
hold off
xlabel("Time(s)")
ylabel("Drift Velocity (m/s)")

figure()
hold on 
grid on
grid minor
box on
plot(time,batt_voltage,'-*','LineWidth',2)
hold off

xlabel("Time(s)")
ylabel("Speed (m/s)")

a = atmosisa(gps_alt/10^3);

figure()
hold on 
grid on
grid minor
box on
plot(time,speed./a,'-*','LineWidth',2)
hold off

xlabel("Time(s)")
ylabel("Mach Number")

figure()
hold on 
grid on
grid minor
box on
plot(time,gps_sat,'-*','LineWidth',2)
hold off

xlabel("Time(s)")
ylabel("Satellites")

figure()
hold on 
grid on
grid minor
box on
plot(time,gps_lat,'-*','LineWidth',2)
hold off

xlabel("Time(s)")
ylabel("Latitude")

m = 0.905; % cansat mass (kg)

fx = m*ax;
fy = m*ay;
fz = m*az;

