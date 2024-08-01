% CLARK Y Parafoil Aerodynamic Analysis
% started 22/07/24 - Rosalind Aves
clear all; close all; clc;

a0 = 6.89 ; %* pi/180; % Clark Y lift curve slope (deg^-1)
alpha_zl = -7 * pi/180; % zero lift angle of attack (deg)

AR(1) = 3;
AR(2) = 3.33;

tau(1) = 0.097; % AR = 3.
tau(2) = 0.105; % AR = 3.33

delta(1) = 0.01875; % AR = 3
delta(2) = 0.022;   % AR = 3.33

a0_dash = (2*pi.*AR) .* tanh(a0./(2*pi.*AR));
a = (pi .* AR .* a0_dash) ./ (pi .* AR + a0_dash .* (1+tau)); % * pi/180; % per degree

alpha = [-7:0.05:10] .* pi/180;
for i = 1:length(a)
    C_L(:,i) = a(i) .* (alpha - alpha_zl);
end

figure();
plot(alpha.*180/pi, C_L(:,1));
hold on
plot(alpha.*180/pi, C_L(:,2));
xlabel('alpha (deg)'); ylabel('C_L');
legend('AR = 3', 'AR = 3.3');
hold off


C_D0 = 0.015 + 0.004 + 0.5 * 0.11 + 0.0001;
for i = 1:length(delta)
    C_D(:,i) = C_D0 + (C_L(:,i).^2 .* (1 + delta(i))) ./ (pi * AR(i));
end

figure();
plot(alpha.*180/pi, C_D(:,1));
hold on
plot(alpha.*180/pi, C_D(:,2));
xlabel('alpha (deg)'); ylabel('C_D');
legend('AR = 3', 'AR = 3.3');
hold off

figure();
plot(alpha.*180/pi, C_L(:,1)./C_D(:,1));
hold on
plot(alpha.*180/pi, C_L(:,2)./C_D(:,2));
xlabel('alpha (deg)'); ylabel('C_L/C_D');
legend('AR = 3', 'AR = 3.3');
hold off

b = 1.42; % span (m)
R_b = linspace(0.4, 1.6, length(alpha));
R = b .* R_b;
beta = b ./ (4.*R);
alpha = 5 * pi/180;

for i = 1:length(a)
    C_L2(:,i) = a(i) .* (alpha - alpha_zl) .* cos(beta).^2;
end

figure();
plot(R_b, C_L2(:,1)./C_L(:,1));
hold on
plot(R_b, C_L2(:,2)./C_L(:,2));
xlabel('line length / span (R/b)'); ylabel('C_L/C_L (R/b=inf)');
legend('AR = 3', 'AR = 3.3');
hold off

S = b^2 ./ AR;
c =  b./AR;
d = 1.5e-3; % line diameter (m)
n = 12;
A_cube = 0.1 * 0.1; % cubesat area (m^2)

for i = 1:length(a)
    C_Dl(:,i) = (n .* R .* d .* cos(alpha).^3) ./ S(i);
    C_Ds(:,i) = (A_cube) ./ S(i);
end

for i = 1:length(delta)
    C_D2(:,i) = C_D0 + C_Dl(:,i)  + C_Ds(:,i)+ (C_L(:,i).^2 .* (1 + delta(i))) ./ (pi * AR(i));
end

figure();
plot(R_b, C_L2(:,1)./C_D2(:,1));
hold on
plot(R_b, C_L2(:,2)./C_D2(:,2));
xlabel('line length / span (R/b)'); ylabel('C_L/C_D');
legend('AR = 3', 'AR = 3.3');
hold off
%% design points

R_b_dp = 0.6; % optimum based on last figure
alpha_dp = 5*pi/180; % before stall
R_dp = R_b_dp * b;
beta_dp = b / (4*R_dp);

C_L2_dp = a(2) * (alpha_dp - alpha_zl) * cos(beta_dp)^2;
C_L_dp = a(2) * (alpha_dp - alpha_zl);

C_Dl_dp = (n * R_dp * d * cos(alpha_dp)^3) / S(2);
C_Ds_dp = (A_cube) / S(2);

C_D2_dp = C_D0 + C_Dl_dp  + C_Ds_dp + (C_L_dp^2 * (1 + delta(2))) / (pi * AR(2));

gamma = atan(1/(C_L2_dp / C_D2_dp)); % glide angle
rho = 1.225;
Wt = 3*9.81;
Dr = (C_L2_dp^2 + C_D2_dp^2)^(0.5);
V = (2*Wt/(rho*S(2)*Dr))^(0.5);
u = V*cos(gamma);
w = V*sin(gamma);
%% 
% Rigging angle stuff

m_s = 3;
g = 9.81;
C_Mc4 = 0;
%mu = -3.35*pi/180;
mu1 = [-10:1:15] .* pi/180;
alpha1 = [-7:0.05:30] .* pi/180;
for j = 1: length(mu1)
    for i = 1 : length(alpha1)
        C_M(j,i) = C_Mc4 - R_dp*(C_Ds_dp*cos(alpha1(1,i) + mu1(1,j))) - R_dp*n*R_dp*d*((cos(alpha1(1,i) + mu1(1,j)))^2)./(S(2)*2*c(2))- m_s*g*R_dp*sin(alpha1(1,i) + mu1(1,j) - gamma)./(0.5*rho*V*V*S(2)*c(2));
    end
end
%%
y=zeros(size(alpha1));
x=alpha_dp.*ones(size(alpha1));
figure();
plot(alpha1,y);
hold on
plot(x,alpha1);
hold on
for j = 1:length(mu1)
    plot(alpha1,C_M(j,:))
end