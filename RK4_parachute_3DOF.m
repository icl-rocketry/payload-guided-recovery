% 3DOF sim
clear; close all
T = 0.05;
end_time = 1000;
num_of_steps = floor(end_time / T);

%% SETUP PARAMETERS
g = 9.81; % m/s

[pfoilParams.b, pfoilParams.c, pfoilParams.S, pfoilParams.AR, pfoilParams.t, pfoilParams.mu, pfoilParams.eps, pfoilParams.a, pfoilParams.R, pfoilParams.d, pfoilParams.n, pfoilParams.m_s, pfoilParams.m_p, pfoilParams.A_cube, ~, pfoilParams.l_cont] = calcPfoilGeometry();
% for i = 1:10
run NoControl.m  % test with no line actuation
% run OneControl.m % test with actuating one control line
% run TwoControl.m % test with actuating both control lines symmetrically
% dxL = -0.0316; % Let out left line 1cm
% dxR = 0.0316;

deltaR = atan(2*dxR / pfoilParams.b); % right line deflection angle
deltaL = atan(2*dxL / pfoilParams.b); % left line deflection angle
deltaA = deltaR - deltaL; % asymmetric deflection angle
deltaS = (deltaR + deltaL)/2; % symmetric deflection angle

% end
u = [deltaS, deltaA]; % control input

[aeroParams] = calcAeroCoeffs(pfoilParams, u);

A = struct2cell(pfoilParams);
pfoilList = cat(2,A{:});
A = struct2cell(aeroParams);
aeroList = cat(2,A{:});
%%

%initial conditions
vel0 = [7.72; 7.87; 50]; % initial ground speed vector
NED0 = [1100; 727; -2625];
W0 = [6.1518; 6.1518; 0]; % initial wind 6.1518
Psi = 0; 
phi = 0; % no initial control
gamma0 = atan(-vel0(3) / vel0(1)); % flight path angle

R_WN = [cos(Psi)*cos(gamma0) sin(Psi)*cos(gamma0) -sin(gamma0); ...
        cos(Psi)*sin(gamma0)*sin(phi) - sin(Psi)*cos(phi) sin(Psi)*sin(gamma0)*sin(phi) + cos(Psi)*cos(phi) cos(gamma0)*sin(phi); ...
        cos(Psi)*sin(gamma0)*cos(phi) + sin(Psi)*sin(phi) sin(Psi)*sin(gamma0)*cos(phi) - cos(Psi)*sin(phi) cos(gamma0)*cos(phi)];
R_BW = eye(3);
Va0 = vel0 - R_BW* R_WN * W0;
Va0 = norm(Va0);

%%

opts = odeset('Events',@iHitTheGround, 'OutputFcn',@odeplot);
[t,x] = ode15s(@(t,x) three_dof_parachute(x, u, W0, aeroParams, pfoilParams, g), [0 end_time], [Va0; gamma0; Psi; 0; 0; NED0(3)], opts);

function [value, isterminal, direction] = iHitTheGround(t,x)
value = x(6);
isterminal = 1;
direction = 0;
end

function [value, isterminal, direction] = cutScript(t,x)
value = x(6);
isterminal = 1;
direction = 0;
end

%% Plots

X = x(:,4);
Y = x(:,5);
Z = -x(:,6);

phi = 0;
theta = x(:,2);
psi = x(:,3);
Va = x(:,1);

for i = 1:size(x,1)
R_BN = [cos(psi(i))*cos(theta(i)) sin(psi(i))*cos(theta(i)) -sin(theta(i)); ...
        cos(psi(i))*sin(theta(i))*sin(phi)-sin(psi(i))*cos(phi) sin(psi(i))*sin(theta(i))*sin(phi)+cos(psi(i))*cos(phi) cos(theta(i))*sin(phi); ...
        cos(psi(i))*sin(theta(i))*cos(phi)+sin(psi(i))*sin(phi) sin(psi(i))*sin(theta(i))*cos(phi)-cos(psi(i))*sin(phi) cos(theta(i))*cos(phi)];

velocityGS(i,:) = eye(3) * [Va(i); 0; 0] + R_BN * W0;

DCM= [1 0 0; 0 cos(phi) sin(phi); ...
    0 -sin(phi) cos(phi)] * ...
    [cos(theta(i)) 0 -sin(theta(i)); 0 1 0; ...
    sin(theta(i)) 0 cos(theta(i))] * ...
    [cos(psi(i)) sin(psi(i)) 0; ...
    -sin(psi(i)) cos(psi(i)) 0; 0 0 1];

velocityI(i,:) = DCM' * velocityGS(i,:)';

end
U = velocityI(:,1);
V = velocityI(:,2);
W = velocityI(:,3);
% V = x(:,8);
% W = x(:,9);

figure();
subplot(3,1,1)
plot(t, X);
hold on;
plot(t, Y);
plot(t, Z);
legend('x', 'y', 'z')
xlabel('time'); ylabel('Distance');
hold off

subplot(3,1,2)
plot(t, U);
hold on
plot(t, V);
plot(t, W);
legend('U', 'V', 'W')
xlabel('time'); ylabel('Velocity');
hold off
% figure();
% plot(t, Cl);
% title("CL vs time");
% xlabel('time'); ylabel('CL');

subplot(3,1,3)
plot(t, sqrt(U.^2 + V.^2 + W.^2));
title("Vmag vs time");
xlabel('time'); ylabel('Vmag');

figure(); 
plot(t, (x(:,2)).*180/pi);
title("Flight Path Angle vs time");
xlabel('time (s)'); ylabel('Angle (deg)');

figure; 
plot3(X,Y,Z); 
axis equal
title("NED");
xlabel('N'); ylabel('E'); zlabel('D');

%% Drift

max_drift_pay = mean([max(X) max(Y)]);
rocket_drift = 1000;

drift = rocket_drift + max_drift_pay;

% U_pretie = U; W_pretie = W; X_pretie = X; Y_pretie = Y; Z_pretie = Z; drift_pretie = drift; save('drift_pretie.mat', 'X_pretie', 'Y_pretie', 'Z_pretie', 'drift_pretie', 'U_pretie', 'W_pretie');