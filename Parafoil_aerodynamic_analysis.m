% CLARK Y Parafoil Aerodynamic Analysis
% started 22/07/24 - Rosalind Aves
clear all; close all; clc;

% Constants and parameters
a0 = rad2deg(0.11) ; % Clark Y 11.7% lift curve slope (deg^-1)
alpha_zl = -3 * pi/180; % Zero lift angle of attack (rad)

% Aspect Ratios
AR = [3, 3.42];

% Tau and delta values for different AR
tau = [0.097, 0.108];
delta = [0.01875, 0.023];

% Calculating a0_dash and a for different AR values
a0_dash = (2 * pi .* AR) .* tanh(a0 ./ (2 * pi .* AR));
a = (pi .* AR .* a0_dash) ./ ((pi .* AR + a0_dash .* (1 + tau)));

% Angle of attack range
alpha = [-3:0.05:10] .* pi/180;

% Lift coefficient calculation
for i = 1:length(a)
    C_L(:,i) = a(i) .* (alpha - alpha_zl);
    C_Ldes(i) = a(i) .* (5*pi/180 - alpha_zl);
end

% Plotting C_L vs alpha
figure();
plot(alpha .* 180/pi, C_L(:,1));
hold on
plot(alpha .* 180/pi, C_L(:,2));
xlabel('alpha (deg)');
ylabel('C_L');
legend('AR = 3', 'AR = 3.42');
hold off

% Drag coefficient calculation
C_D0 = 0.015 + 0.004 + 0.5 * 0.117 + 0.0001;
for i = 1:length(delta)
    C_D(:,i) = C_D0 + (C_L(:,i).^2 .* (1 + delta(i))) ./ (pi * AR(i));
    C_Ddes(i) = C_D0 + (C_Ldes(i).^2 .* (1 + delta(i))) ./ (pi * AR(i));
end

% Plotting C_D vs alpha
figure();
plot(alpha .* 180/pi, C_D(:,1));
hold on
plot(alpha .* 180/pi, C_D(:,2));
xlabel('alpha (deg)');
ylabel('C_D');
legend('AR = 3', 'AR = 3.42');
hold off

% Plotting C_L/C_D vs alpha
figure();
plot(alpha .* 180/pi, C_L(:,1) ./ C_D(:,1));
hold on
plot(alpha .* 180/pi, C_L(:,2) ./ C_D(:,2));
xlabel('alpha (deg)');
ylabel('C_L/C_D');
legend('AR = 3', 'AR = 3.42');
hold off

% Span and line length ratio calculations
b = 1.27; % Span (m)
R_b = linspace(0.4, 1.6, 100);
R = b .* R_b;
R_des = 1.2;
beta = b ./ (4 .* R);
beta_des = b./(4*R_des);
alpha_des = 5 * pi/180;

for i = 1:length(a)
    C_L2(:,i) = a(i) .* (alpha_des - alpha_zl) .* cos(beta).^2;
    C_L3(:,i) = a(i) .* (alpha - alpha_zl) .* cos(beta_des).^2;
    C_L2des(:,i) = a(i) .* (alpha_des - alpha_zl) .* cos(beta_des).^2;
end

% Plotting C_L2/C_L vs R/b
figure();
plot(R_b, C_L2(:,1) ./ C_Ldes(1));
hold on
plot(R_b, C_L2(:,2) ./ C_Ldes(2));
xlabel('line length / span (R/b)');
ylabel('C_L/C_L (R/b=inf)');
legend('AR = 3', 'AR = 3.3');
hold off

% Surface area and chord length calculations
S = b^2 ./ AR;
c = b ./ AR;
d = 1.5e-3; % Line diameter (m)
n = 24;
A_cube = 0.1 * 0.1; % Cubesat area (m^2)

for i = 1:length(a)
    C_Dl_des(:,i) = (n .* R .* d .* cos(alpha_des).^3) ./ S(i);
    C_Ds(:,i) = (A_cube) ./ S(i);
end

for i = 1:length(delta)
    C_D2(:,i) = C_D0 + C_Dl_des(:,i) + C_Ds(:,i) + (C_Ldes(:,i).^2 .* (1 + delta(i))) ./ (pi * AR(i));
end

% Plotting C_L2/C_D2 vs R/b
figure();
plot(R_b, C_L2(:,1) ./ C_D2(:,1));
hold on
plot(R_b, C_L2(:,2) ./ C_D2(:,2));
xlabel('line length / span (R/b)');
ylabel('C_L/C_D');
legend('AR = 3', 'AR = 3.42');
hold off

%% Design points

R_b_dp = 1.2; % Based on line lengths
alpha_dp = 2 * pi/180; % Before stall
R_dp = R_b_dp * b;
beta_dp = b / (4 * R_dp);

C_L2_dp = a(2) * (alpha_dp - alpha_zl) * cos(beta_dp)^2;
C_L_dp = a(2) * (alpha_dp - alpha_zl);

C_Dl_dp = (n * R_dp * d * cos(alpha_dp)^3) / S(2);
C_Ds_dp = (A_cube) / S(2);

C_D2_dp = C_D0 + C_Dl_dp + C_Ds_dp + (C_L_dp^2 * (1 + delta(2))) / (pi * AR(2));

gamma = atan(1 / (C_L2_dp / C_D2_dp)); % Glide angle
rho = 1.225; % Air density (kg/m^3)
Wt = (3+0.1) * 9.81; % Weight (N)
Ct = sqrt(C_L2_dp^2 + C_D2_dp^2);
V = sqrt((2 * Wt) / (rho * S(2) * Ct));
u = V * cos(gamma);
w = V * sin(gamma);

%% Rigging angle analysis

m_s = 3; % Mass (kg)
g = 9.81; % Gravity (m/s^2)
C_Mc4 = -0.08; %Cm0
mu1 = [0:1:16] .* pi/180;
alpha1 = [-3:0.05:30] .* pi/180;

for j = 1:length(mu1)
    for i = 1:length(alpha1)
        C_M(j,i) = C_Mc4 - R_dp/c(2) * (C_Ds_dp * cos(alpha1(i) + mu1(j))) - R_dp * n * R_dp * d * (cos(alpha1(i) + mu1(j))^2) / (S(2) * 2 * c(2)) - m_s * g * R_dp * sin(alpha1(i) + mu1(j) - gamma) / (0.5 * rho * V^2 * S(2) * c(2));
    end
end

% Plotting rigging angle analysis
y = zeros(size(alpha1));
x = alpha_dp * ones(size(alpha1));
figure();
plot(rad2deg(alpha1), y);
hold on
plot(x, (alpha1));
hold on
for j = 1:length(mu1)
    plot(rad2deg(alpha1), C_M(j,:));
end
xlabel('alpha (deg)');
ylabel('C_M');
legend('Zero Line', 'Design Point', 'C_M curves');
hold off


%% Lateral dynamics

C_nr = -C_D2_dp/6;
C_ncl = 0.01;

r = -C_ncl/C_nr * (2*V)/b * delta(2);
R_c = u/r;

f = asind((u*r)/g);

%% curve fitting aero coeffs

%y = C_D
% C_D(:,i) = C_D0 + (C_L(:,i).^2 .* (1 + delta(i))) ./ (pi * AR(i));
%x = alpha
alpha = [-3:0.05:10] .* pi/180;

pCD = fit(alpha', C_D(:,2),'poly1');
y1 = pCD(alpha);

%y = C_L
%x = alpha
pCL = fit(alpha', C_L(:,2),'poly1');
% y2 = polyval(pCL, alpha);

%y = C_M
%x = alpha
alpha1 = [-3:0.05:30] .* pi/180;
pCM = fit(alpha1', C_M(end,:)','poly1');
% y3 = polyval(pCM, alpha1);

% figure
% plot(alpha,C_D(:,2),'o')
% hold on
% plot(alpha,y1)
% hold off