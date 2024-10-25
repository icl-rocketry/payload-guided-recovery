clear; close all
T = 0.05;
end_time = 1000;
num_of_steps = floor(end_time / T);

%% SETUP PARAMETERS
g = 9.81; % m/s

[pfoilParams.b, pfoilParams.c, pfoilParams.S, pfoilParams.AR, pfoilParams.t, pfoilParams.mu, pfoilParams.eps, pfoilParams.a, pfoilParams.R, pfoilParams.d, pfoilParams.n, pfoilParams.m_s, pfoilParams.m_p, pfoilParams.A_cube, ~, pfoilParams.l_cont] = calcPfoilGeometry();

run NoControl.m  % test with no line actuation
% run OneControl.m % test with actuating one control line
% run TwoControl.m % test with actuating both control lines symmetrically
% pfoilParams.eps = deg2rad(20);
deltaR = atan(2*dxR / pfoilParams.b); % right line deflection angle
deltaL = atan(2*dxL / pfoilParams.b); % left line deflection angle
deltaA = deltaR - deltaL; % asymmetric deflection angle
deltaS = (deltaR + deltaL)/2; % symmetric deflection angle

u = [deltaS deltaA]; % control input

[aeroParams] = calcAeroCoeffs(pfoilParams, u);

% aeroParams.CLalpha = 2.56;

%%

opts = odeset('Events',@iHitTheGround, 'OutputFcn',@odeplot);
[t,x] = ode23s(@(t,x) six_dof_parachute(x, u, aeroParams, pfoilParams, g), [0 end_time], [0; 0; 3000; 0; 0; 0; 0; 0; 45; 0; 0; 0], opts);

function [value, isterminal, direction] = iHitTheGround(t,x)
value = x(3);
isterminal = 1;
direction = 0;
end

% for i = 1:num_of_steps
%     % x{step}(1:3)=[0;0;0];
%     x{i+1}=RungeKutta4(@six_dof_parachute, x{i}, u{i}, T);
%
%     if disp == 100
%         disp = 0;
%         T*i;
%     else
%         disp = disp + 1;
%     end
%
%
%     % x(1)=X (North-East-Up format)
%     % x(2)=Y
%     % x(3)=Z
%     % x(4)=roll (phi)
%     % x(5)=pitch (theta)
%     % x(6)=yaw (psi)
%     % x(7)=d/dt X
%     % x(8)=d/dt Y
%     % x(9)=d/dt Z
%     % x(10)=d/dt roll
%     % x(11)=d/dt pitch
%     % x(12)=d/dt yaw
%     %
%     % u(1)=brake
%     % u(2)=aileron
%     X(i) = x{i}(1);
%     Y(i) = x{i}(2);
%     Z(i) = x{i}(3);
%
%     phi(i) = x{i}(4);
%     theta(i) = x{i}(5);
%     psi(i) = x{i}(6);
%
%     U(i) = x{i}(7);
%     V(i) = x{i}(8);
%     W(i) = x{i}(9);
%
%     P(i) = x{i}(10);
%     Q(i) = x{i}(11);
%     R(i) = x{i}(12);
%
%     t(i) = i .* T;
%
%     if x{1}(3) ~= 0
%         if x{i+1}(3) <= 0
%             break;
%         end
%     end
% end
%
%
%
%
% figure();
% plot(t, X);
% title("X vs time");
% xlabel('time'); ylabel('X');
%
% figure();
% plot(t, U);
% title("U vs time");
% xlabel('time'); ylabel('U');
%
% figure();
% plot(t, Z);
% title("Z vs time");
% xlabel('time'); ylabel('Z');
% % figure();
% % plot3(X, Y, Z);
%
% figure();
% plot(t, W);
% title("W vs time");
% xlabel('time'); ylabel('W');
%
% % figure();
% % plot(t, Cl);
% % title("CL vs time");
% % xlabel('time'); ylabel('CL');


%%
X = x(:,1);
Y = x(:,2);
Z = x(:,3);

phi = x(:,4);
theta = x(:,5);
psi = x(:,6);

for i = 1:size(x,1)
DCM= [1 0 0; 0 cos(phi(i)) sin(phi(i)); ...
    0 -sin(phi(i)) cos(phi(i))] * ...
    [cos(theta(i)) 0 -sin(theta(i)); 0 1 0; ...
    sin(theta(i)) 0 cos(theta(i))] * ...
    [cos(psi(i)) sin(psi(i)) 0; ...
    -sin(psi(i)) cos(psi(i)) 0; 0 0 1];



velocityI(i,1:3) = DCM' * [x(i,7); x(i,8); x(i,9)];


end
U = velocityI(:,1);
V = velocityI(:,2);
W = velocityI(:,3);
% V = x(:,8);
% W = x(:,9);

P = x(:,10);
Q = x(:,11);
R = x(:,12);

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


% figure(); plot(t, x(:, 10)); hold on; plot(t, x(:, 11)); plot(t, x(:, 12));


%% Make full drift plots

