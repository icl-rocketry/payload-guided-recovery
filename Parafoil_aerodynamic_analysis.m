% CLARK Y Parafoil Aerodynamic Analysis
% started 22/07/24 - Rosalind Aves
clear all; close all; clc;

% Constants and parameters
a0 = 6.89 ; % Clark Y lift curve slope (deg^-1)
alpha_zl = -7 * pi/180; % Zero lift angle of attack (rad)

% Aspect Ratios
AR = [3, 3.33];

% Tau and delta values for different AR
tau = [0.097, 0.105];
delta = [0.01875, 0.022];

% Calculating a0_dash and a for different AR values
a0_dash = (2 * pi .* AR) .* tanh(a0 ./ (2 * pi .* AR));
a = (pi .* AR .* a0_dash) ./ (pi .* AR + a0_dash .* (1 + tau));

% Angle of attack range
alpha = [-7:0.05:10] .* pi/180;

% Lift coefficient calculation
for i = 1:length(a)
    C_L(:,i) = a(i) .* (alpha - alpha_zl);
end

% Plotting C_L vs alpha
figure();
plot(alpha .* 180/pi, C_L(:,1));
hold on
plot(alpha .* 180/pi, C_L(:,2));
xlabel('alpha (deg)');
ylabel('C_L');
legend('AR = 3', 'AR = 3.3');
hold off

% Drag coefficient calculation
C_D0 = 0.015 + 0.004 + 0.5 * 0.11 + 0.0001;
for i = 1:length(delta)
    C_D(:,i) = C_D0 + (C_L(:,i).^2 .* (1 + delta(i))) ./ (pi * AR(i));
end

% Plotting C_D vs alpha
figure();
plot(alpha .* 180/pi, C_D(:,1));
hold on
plot(alpha .* 180/pi, C_D(:,2));
xlabel('alpha (deg)');
ylabel('C_D');
legend('AR = 3', 'AR = 3.3');
hold off

% Plotting C_L/C_D vs alpha
figure();
plot(alpha .* 180/pi, C_L(:,1) ./ C_D(:,1));
hold on
plot(alpha .* 180/pi, C_L(:,2) ./ C_D(:,2));
xlabel('alpha (deg)');
ylabel('C_L/C_D');
legend('AR = 3', 'AR = 3.3');
hold off

% Span and line length ratio calculations
b = 1.42; % Span (m)
R_b = linspace(0.4, 1.6, length(alpha));
R = b .* R_b;
beta = b ./ (4 .* R);
alpha = 5 * pi/180;

for i = 1:length(a)
    C_L2(:,i) = a(i) .* (alpha - alpha_zl) .* cos(beta).^2;
end

% Plotting C_L2/C_L vs R/b
figure();
plot(R_b, C_L2(:,1) ./ C_L(:,1));
hold on
plot(R_b, C_L2(:,2) ./ C_L(:,2));
xlabel('line length / span (R/b)');
ylabel('C_L/C_L (R/b=inf)');
legend('AR = 3', 'AR = 3.3');
hold off

% Surface area and chord length calculations
S = b^2 ./ AR;
c = b ./ AR;
d = 1.5e-3; % Line diameter (m)
n = 12;
A_cube = 0.1 * 0.1; % Cubesat area (m^2)

for i = 1:length(a)
    C_Dl(:,i) = (n .* R .* d .* cos(alpha).^3) ./ S(i);
    C_Ds(:,i) = (A_cube) ./ S(i);
end

for i = 1:length(delta)
    C_D2(:,i) = C_D0 + C_Dl(:,i) + C_Ds(:,i) + (C_L(:,i).^2 .* (1 + delta(i))) ./ (pi * AR(i));
end

% Plotting C_L2/C_D2 vs R/b
figure();
plot(R_b, C_L2(:,1) ./ C_D2(:,1));
hold on
plot(R_b, C_L2(:,2) ./ C_D2(:,2));
xlabel('line length / span (R/b)');
ylabel('C_L/C_D');
legend('AR = 3', 'AR = 3.3');
hold off

%% Design points

R_b_dp = 0.6; % Optimum based on the last figure
alpha_dp = 5 * pi/180; % Before stall
R_dp = R_b_dp * b;
beta_dp = b / (4 * R_dp);

C_L2_dp = a(2) * (alpha_dp - alpha_zl) * cos(beta_dp)^2;
C_L_dp = a(2) * (alpha_dp - alpha_zl);

C_Dl_dp = (n * R_dp * d * cos(alpha_dp)^3) / S(2);
C_Ds_dp = (A_cube) / S(2);

C_D2_dp = C_D0 + C_Dl_dp + C_Ds_dp + (C_L_dp^2 * (1 + delta(2))) / (pi * AR(2));

gamma = atan(1 / (C_L2_dp / C_D2_dp)); % Glide angle
rho = 1.225; % Air density (kg/m^3)
Wt = 3 * 9.81; % Weight (N)
Dr = sqrt(C_L2_dp^2 + C_D2_dp^2);
V = sqrt((2 * Wt) / (rho * S(2) * Dr));
u = V * cos(gamma);
w = V * sin(gamma);

%% Rigging angle analysis

m_s = 3; % Mass (kg)
g = 9.81; % Gravity (m/s^2)
C_Mc4 = 0;
mu1 = [-10:1:15] .* pi/180;
alpha1 = [-7:0.05:30] .* pi/180;

for j = 1:length(mu1)
    for i = 1:length(alpha1)
        C_M(j,i) = C_Mc4 - R_dp * (C_Ds_dp * cos(alpha1(i) + mu1(j))) - R_dp * n * R_dp * d * (cos(alpha1(i) + mu1(j))^2) / (S(2) * 2 * c(2)) - m_s * g * R_dp * sin(alpha1(i) + mu1(j) - gamma) / (0.5 * rho * V^2 * S(2) * c(2));
    end
end

% Plotting rigging angle analysis
y = zeros(size(alpha1));
x = alpha_dp * ones(size(alpha1));
figure();
plot(alpha1, y);
hold on
plot(x, alpha1);
hold on
for j = 1:length(mu1)
    plot(alpha1, C_M(j,:));
end
xlabel('alpha (deg)');
ylabel('C_M');
legend('Zero Line', 'Design Point', 'C_M curves');
hold off
