% CLARK Y Parafoil Aerodynamic Analysis
% started 22/07/24 - Rosalind Aves
clear all; close all; clc;

% Constants and parameters
a0 = 0.11*180/pi; % Clark Y lift curve slope (deg^-1)
alpha_zl = -3 * pi/180; % Zero lift angle of attack (rad)

% Aspect Ratios
AR = [3, 3.42];

% Tau and delta values for different AR
tau = [0.097, 0.108];
delta = [0.01875, 0.023];

% Calculating a0_dash and a for different AR values
a0_dash = (2 * pi .* AR) .* tanh(a0 ./ (2 * pi .* AR));
a = (pi .* AR .* a0_dash) ./ (pi .* AR + a0_dash .* (1 + tau));
% Angle of attack range
alpha = [-3:0.05:10] .* pi/180;% rad

% Lift coefficient calculation
for i = 1:length(a)
    C_L(:,i) = a(i) .* (alpha - alpha_zl);
end
% Lift coefficient calculation
for i = 1:length(a)
    C_Lalpha(:,i) = a(i);
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
%% 

% Span and line length ratio calculations
b = 1.27; % Span (m)
R_b = linspace(0.4, 1.6, length(alpha));
R = b .* R_b;
beta = b ./ (4 .* R);
% alpha = 5 * pi/180;
e0 = 30*pi/180; % rad
e = 0.5;
for i = 1:length(a)
    C_L2(:,i) = a(i) .* (alpha - alpha_zl) .* cos(beta).^2;
    k(:,i) = AR(i)*pi/a0;
    k1(:,i) = (sqrt(k(:,i)^2 + 1) + 1)/(sqrt(k(:,i)^2 + 4) + 2);
    k2(:,i) = (sqrt(k(:,i)^2 + 4) - 1)/(sqrt(k(:,i)^2 + 1) + 1);
    C_L3alpha(:,i) = AR(i)*pi/(sqrt(k(i)^2 + 1)+1); 
    C_L3(:,i) = C_L3alpha(:,i)*(alpha*cos(e0/2)-alpha_zl)*cos(e0/2);
end
% Plotting C_L2 vs alpha
figure();
xline(0);
hold on
plot(alpha .* 180/pi, C_L2(:,1));
hold on
plot(alpha .* 180/pi, C_L2(:,2));
hold on
plot(alpha .* 180/pi, C_L3(:,1));
hold on
plot(alpha .* 180/pi, C_L3(:,2));
xlabel('alpha');
ylabel('C_{L2} & C_{L3}');
legend('alpha = 0','AR = 3', 'AR = 3.42', 'C_{L3} AR = 3 ','C_{L3} AR = 3.42 ');
hold off

for i = 1:length(a)
    C_L2alpha(:,i) = a(i).* cos(beta).^2;%per rad
end
% Plotting C_L2alpha vs alpha
figure();
plot(alpha .* 180/pi, C_L2alpha(:,1));
hold on
plot(alpha .* 180/pi, C_L2alpha(:,2));
hold on
plot(alpha .* 180/pi, C_L3alpha(:,1)*ones(size(alpha)));
hold on
plot(alpha .* 180/pi, C_L3alpha(:,2)*ones(size(alpha)));
hold on
xline(0);
xlabel('alpha');
ylabel('C_{L2_{alpha}} & C_{L3_{alpha}}');
legend('C_{L2_{alpha}} AR = 3', 'C_{L2_{alpha}} AR = 3.42','C_{L3_{alpha}} AR = 3 ','C_{L3_{alpha}} AR = 3.42 ', 'alpha = 0');
hold off

% Plotting C_L2/C_L and C_L3/C_L vs R/b
figure();
plot(R_b, C_L2(:,1) ./ C_L(:,1));
hold on
plot(R_b, C_L2(:,2) ./ C_L(:,2));
hold on
plot(R_b, C_L3(:,1) ./ C_L(:,1));
hold on
plot(R_b, C_L3(:,2) ./ C_L(:,2));
xlabel('line length / span (R/b)');
ylabel('C_L/C_L (R/b=inf)');
legend('C_{L2}/C_L AR = 3', 'C_{L2}/C_L AR = 3.42', 'C_{L3}/C_L AR = 3', 'C_{L3}/C_L AR = 3.42');
hold off

% Surface area and chord length calculations
S = b^2 ./ AR;
c = b ./ AR;
d = 1.5e-3; % Line diameter (m)
n = 24;
A_cube = sqrt(0.3^2 + 0.1^2)*0.1; % Cubesat area (m^2)

for i = 1:length(a)
    C_Dl(:,i) = (n .* R .* d .* cos(alpha).^3) ./ S(i);
    C_Ds(:,i) = (A_cube) ./ S(i);
end

for i = 1:length(delta)
    C_D2(:,i) = C_D0 + C_Dl(:,i) + C_Ds(:,i) + (C_L(:,i).^2 .* (1 + delta(i))) ./ (pi * AR(i));
    C_D3(:,i) = C_D0 + + C_Ds(:,i) + C_L3alpha(:,i).^2 * (alpha*cos(e0/2) - alpha_zl).^2 / (e*AR(:,i)*pi);
end

for i = 1:length(a)
    C_Dlalpha(:,i) = -(n .* R .* d .*3.* (alpha).^2.*sin(alpha).^3) ./ S(i);
    % C_Dlalpha(:,i) = C_Dlalpha(:,i) * pi/180;
end

for i = 1:length(delta)
    C_D2alpha(:,i) =  C_Dlalpha(:,i) + (2.*C_L(:,i).*C_Lalpha(:,i) .* (1 + delta(i))) ./ (pi * AR(i));
    % C_D2alpha(:,i) = C_D2alpha(:,i) * pi/180;
    C_D3alpha(:,i) = C_L3alpha(:,i)^2 * 2 * cos(e0/2)* (alpha*cos(e0/2) - alpha_zl)/ (e*AR(:,i)*pi);
end
% Plotting C_D2 and C_D3 vs alpha
figure();
xline(0);
hold on
plot(alpha .* 180/pi, C_D3(:,1));
hold on
plot(alpha .* 180/pi, C_D3(:,2));
hold on
plot(alpha .* 180/pi, C_D2(:,1));
hold on
plot(alpha .* 180/pi, C_D2(:,2));
xlabel('alpha');
ylabel('C_{D2} and C_{D3}');
legend('alpha = 0','C_{D3} AR = 3', 'C_{D3} AR = 3.42','C_{D2} AR = 3', 'C_{D2} AR = 3.42');
hold off

% Plotting C_D2alpha and C_D3alpha vs alpha
figure();
xline(0);
hold on
plot(alpha .* 180/pi, C_D3alpha(:,1));
hold on
plot(alpha .* 180/pi, C_D3alpha(:,1));
hold on
plot(alpha .* 180/pi, C_D2alpha(:,1));
hold on
plot(alpha .* 180/pi, C_D2alpha(:,2));
xlabel('alpha');
ylabel('C_{D2_{alpha}} and C_{D3_{alpha}}');
legend('alpha = 0','C_{D3_{alpha}} AR = 3', 'C_{D3_{alpha}} AR = 3.42','C_{D2_{alpha}} AR = 3', 'C_{D2_{alpha}} AR = 3.42');
hold off

% Plotting C_L2/C_D2 vs R/b
figure();
plot(alpha*180/pi, C_L2(:,1) ./ C_D2(:,1));
hold on
plot(alpha*180/pi, C_L2(:,2) ./ C_D2(:,2));
xlabel('alpha(deg)');
ylabel('C_L/C_D');
legend('AR = 3', 'AR = 3.42');
hold off

% Plotting C_L2/C_D2 vs R/b
figure();
plot(R_b, C_L2(:,1) ./ C_D2(:,1));
hold on
plot(R_b, C_L2(:,2) ./ C_D2(:,2));
xlabel('line length / span (R/b)');
ylabel('C_L/C_D');
legend('AR = 3', 'AR = 3.42');
hold off


%% 

% C_M calculations
% maintaining stable equilibrium at each alpha
rho = 1.225; % Air density (kg/m^3)
R_b_dp = 1.2; % Optimum based on the last figure
R_dp = R_b_dp * b;
m_s = 3 + 0.1; 
g = 9.81; 
C_Mc4 = -0.08;
mu = [-100:0.05:100] .* pi/180; 
corr_mu = NaN(1, length(alpha)); 

for i = 1:length(alpha)
    C_M3q(:,i) = - C_L3alpha(2) * (cos(e0/2)^2)/12;%pitch damping derivative assuming that the CP and reference point coincide
    gamma(:,i) = atan(1/(C_L2(i,2)/C_D2(i,2)));
    Wt = 3 * 9.81; 
    Dr = sqrt(C_L2(i,2)^2 + C_D2(i,2)^2);
    V = sqrt((2 * Wt) / (rho * S(2) * Dr));
    LbyD_alpha(i) = (C_L2alpha(i,2) * C_D2(i,2) - C_L2(i,2) * C_D2alpha(i,2)) / (C_D2(i,2)^2);
    d_gamma_alpha = (1 / (1 + (1 / (C_L2(i,2) / C_D2(i,2)))^2)) * LbyD_alpha(i);
    for j = 1:length(mu)
        C_M2(j,i) = C_Mc4 - (R_dp/c(2)) * (C_Ds(2) * cos(alpha(i) + mu(j))) ...
            - R_dp * n * R_dp * d * (cos(alpha(i) + mu(j))^2) / (S(2) * 2 * c(2)) ...
            - m_s * g * R_dp * sin(alpha(i) + mu(j) - gamma(:,i)) / (0.5 * rho * V^2 * S(2) * c(2)); 

        [~, idx] = min(abs(C_M2(:,i))); 
        corr_mu(i) = mu(idx) * 180/pi; 
        corr_mu(i) = corr_mu(i) * pi/180;
        C_M2_alpha(i) = (R_dp/c(2)) * (C_Ds(2) * sin(alpha(i) + corr_mu(i))) ...
            + R_dp * n * R_dp * d * (sin(alpha(i) + corr_mu(i))^2) / (S(2) * 2 * c(2)) ...
            - d_gamma_alpha * m_s * g * R_dp * cos(alpha(i) + corr_mu(i) - gamma(:,i)) / (0.5 * rho * V^2 * S(2) * c(2));
        % C_M2_alpha(i) = C_M2_alpha(i) * pi/180;
        % Penalize if C_M2_alpha is not negative
        if C_M2_alpha >= 0
            corr_mu(i) = Inf; 
        end
    end
end
% Plotting stable equilibrium mu vs alpha
figure();
xline(0);
hold on
plot(alpha .* 180/pi, corr_mu.*180/pi);
xlabel('alpha');
ylabel('rigging angle for stable equilibrium');
hold off
%updating C_M
for i =1:length(alpha)
    up_C_M2(i) = C_Mc4 - (R_dp/c(2)) * (C_Ds(2) * cos(alpha(i) + corr_mu(i))) ...
            - R_dp * n * R_dp * d * (cos(alpha(i) + corr_mu(i))^2) / (S(2) * 2 * c(2)) ...
            - m_s * g * R_dp * sin(alpha(i) + corr_mu(i) - gamma(:,i)) / (0.5 * rho * V^2 * S(2) * c(2)); 
end
% Plotting updated C_M2 vs alpha
figure();
xline(0);
hold on
plot(alpha .* 180/pi, up_C_M2);
hold on
plot(alpha .* 180/pi, C_M3);
xlabel('alpha');
ylabel('updated C_{M2} and C_{M3}');
legend('alpha = 0','C_{M2}', 'C_{M3}')
hold off

%% data storing

headers = {'alpha (deg)', 'C_L2', 'C_L2alpha', 'C_D2', 'C_D2alpha', 'C_L2/C_D2','up_C_M2', 'C_M2_alpha','corr_mu'};
data = [alpha' .* 180/pi, C_L2(:,2), C_L2alpha(:,2), C_D2(:,2), C_L2(:,2) ./ C_D2(:,2), C_D2alpha(:,2), up_C_M2', C_M2_alpha',corr_mu' .*180/pi]; 
filename = 'aerodynamic_analysis_results.xlsx';
writecell(headers, filename, 'Sheet', 1, 'Range', 'A1');
writematrix(data, filename, 'Sheet', 1, 'Range', 'A2');
%% 
C_nr3 = -(C_D0 * sin(e0))/(3*e0) - (C_L3alpha(2)^2 * sin(e0/2)^2 * alpha_zl)/(e*AR(2)*8*e0) - (C_L3alpha(2)*e0^2)/(AR^2 * 24);
C_lp3 = - C_L3alpha(2) * k1(2) * sin(e0)/(8*e0);
C_lphi3 = 0;
%% Design points

R_b_dp = 1.2; % Optimum based on the last figure
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
C_Mc4 = -0.08;
mu1 = [0:1:16] .* pi/180;
alpha1 = [-3:0.05:30] .* pi/180;

for j = 1:length(mu1)
    for i = 1:length(alpha1)
        C_M(j,i) = C_Mc4 - (R_dp/c(2)) * (C_Ds_dp * cos(alpha1(i) + ...
            mu1(j))) - R_dp * n * R_dp * d * (cos(alpha1(i) + mu1(j))^2) / (S(2) * 2 * c(2)) - m_s * g * R_dp * sin(alpha1(i) + mu1(j) - gamma) / (0.5 * rho * V^2 * S(2) * c(2));
    end
end

% Plotting rigging angle analysis
y = zeros(size(alpha1));
x = (alpha_dp*180/pi) * ones(size(alpha1));
figure();
plot(alpha1.*180/pi, y);
hold on
plot(x, alpha1);
hold on
for j = 1:length(mu1)
    plot(alpha1.*180/pi, C_M(j,:));
end
xlabel('alpha (deg)');
ylabel('C_M')
legend('Zero Line', 'Design Point', 'C_M curves');
hold off
