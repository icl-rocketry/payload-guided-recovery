% CLARK Y Parafoil Aerodynamic Analysis DESIGN POINTS
% NEED TO CONVERT TO DEGREES
% started 25/09/24 - Rosalind Aves

clear all; close all; clc;

% Constants and parameters
a0 = 0.11*180/pi; % Clark Y lift curve slope (deg^-1)
alpha_zl = -3 * pi/180; % Zero lift angle of attack (rad)
alpha = [-3:0.05:10] .* pi/180;% rad
% Aspect Ratios
AR = 3.42;

% Span and line length ratio calculations
b = 1.27; % Span (m)
controlLengthWT = 0.365; % line length for maxRadius (function of)

[eps, Pheight, R] = calcEpsilon(b, controlLengthWT);

beta = b ./ (4 .* R);

% Surface area and chord length calculations
S = b^2 ./ AR;
c = b ./ AR;
d = 1.5e-3; % Line diameter (m)
n = 24;
m_s = 3 + 0.1;
g = 9.81;
A_cube = sqrt(0.3^2 + 0.1^2)*0.1; % Cubesat area (m^2)
rho = 1.225; % Air density (kg/m^3)

%% C_L coefficients !!
tau = 0.108;
a0_dash = (2 * pi .* AR) .* tanh(a0 ./ (2 * pi .* AR));
C_Lalpha_paper1 = (pi .* AR .* a0_dash) ./ (pi .* AR + a0_dash .* (1 + tau)); % FOR AIRFOIL
C_L2alpha_paper1 = C_Lalpha_paper1.* cos(beta).^2; % FOR PARAFOIL per rad

C_L_paper1(:) = C_Lalpha_paper1 .* (alpha - alpha_zl) .* cos(beta).^2; %FOR PARAFOIL per rad
C_L0_paper1 = C_Lalpha_paper1 .* (0 - alpha_zl) .* cos(beta).^2;

k_paper2 = AR*pi/a0;
k1_paper2 = (sqrt(k_paper2^2 + 1) + 1)/(sqrt(k_paper2^2 + 4) + 2);
k2_paper2 = (sqrt(k_paper2^2 + 4) - 1)/(sqrt(k_paper2^2 + 1) + 1);

C_Lalpha_paper2 = AR*pi/(sqrt(k_paper2^2 + 1)+1);
C_L_paper2(:) = C_Lalpha_paper2.*(alpha.*cos(eps/2)-alpha_zl).*cos(eps/2);
C_L0_paper2 = C_Lalpha_paper2.*(0 .*cos(eps/2)-alpha_zl).*cos(eps/2);


%% C_D coefficients !! - no derivatives just pure CD

delta = [0.023]; e = 0.9;
C_D0_paper1 = 0.015 + 0.004 + 0.5 * 0.117 + 0.0001;
C_D_paper1 = C_D0_paper1 + ((C_Lalpha_paper1 * (alpha - alpha_zl)).^2 .* (1 + delta)) ./ (pi * AR);

C_Dl = (n .* R .* d .* cos(alpha).^3) ./ S;
C_Ds = A_cube/ S;

C_D_paper1 = C_D_paper1 + C_Dl + C_Ds; % ADDING LINE AND STORE DRAG
C_D_paper2 = C_D0_paper1 + C_Ds + C_Lalpha_paper2.^2 .* (alpha*cos(eps/2) - alpha_zl).^2 / (e*AR*pi);

%% C_M coefficients !!

% maintaining stable equilibrium at each alpha
C_Mc4 = -0.08; % C_M0
mu = 5 .* pi/180; % imposed incidence by line lengths

C_Mq_paper2 = - C_Lalpha_paper2 .* (cos(eps/2)^2)/12;%pitch damping derivative assuming that the CP and reference point coincide

gamma = atan(1./(C_L_paper1 ./ C_D_paper1));

Wt = m_s * g;
Ct = sqrt(C_L_paper1.^2 + C_D_paper1.^2);
V = sqrt((2 * Wt) ./ (rho .* S .* Ct));
u = V .* cos(gamma);
w = V .* sin(gamma);

C_Ll = - (n * d * R * cos(alpha + mu).^3 .* sin(alpha + mu)) ./ (S);
C_Ls = 0;

C_M_paper1 = C_Mc4 - (R/c) * (C_Ds * cos(alpha + mu)) ...
            - R./ (S * 2 * c) * (n * R * d .* (cos(alpha + mu).^2))  ...
            - (m_s * g * R * sin(alpha + mu - gamma)) ./ (0.5 .* rho * V.^2 .* S .* c);

%% Yaw damping

C_nr_paper2 = -(C_D0_paper1 * sin(eps))/(3*eps) - (C_Lalpha_paper2^2 * sin(eps/2)^2 * alpha_zl)/(e*AR*pi*8*eps^2) - (C_Lalpha_paper2 * eps^2) / (AR^2 * 24);

%% Roll damping

C_lp_paper2 = - C_Lalpha_paper2 * k1_paper2 * sin(eps)/(8*eps);
C_lphi_paper2 = 0; % no value 

%% BINNNED


% for i = 1:length(C_Lalpha_paper1)
%     C_Dlalpha(:,i) = -(n .* R .* d .*3.* (alpha).^2.*sin(alpha).^3) ./ S(i);
%     % C_Dlalpha(:,i) = C_Dlalpha(:,i) * pi/180;
% end
% 
% for i = 1:length(delta)
%     C_D2alpha(:,i) =  C_Dlalpha(:,i) + (2.*C_L(:,i).*C_Lalpha(:,i) .* (1 + delta(i))) ./ (pi * AR(i));
%     % C_D2alpha(:,i) = C_D2alpha(:,i) * pi/180;
%     C_D3alpha(:,i) = C_Lalpha_paper2(:,i)^2 * 2 * cos(eps/2)* (alpha*cos(eps/2) - alpha_zl)/ (e*AR(:,i)*pi);
% end


% C_Mc4 = -0.08;
% mu = [-100:0.05:100] .* pi/180;
% corr_mu = NaN(1, length(alpha));
% 
% for i = 1:length(alpha)
%     C_M3q(:,i) = - C_Lalpha_paper2 * (cos(eps/2)^2)/12;%pitch damping derivative assuming that the CP and reference point coincide
%     gamma(:,i) = atan(1/(C_L_paper1(i)/C_D2(i)));
%     Wt = 3 * 9.81;
%     Ct = sqrt(C_L_paper1(i)^2 + C_D2(i)^2);
%     V = sqrt((2 * Wt) / (rho * S * Ct));
%     L_D_alpha(i) = (C_L2alpha(i) * C_D2(i) - C_L_paper1(i) * C_D2alpha(i)) / (C_D2(i)^2);
%     d_gamma_alpha = (1 / (1 + (1 / (C_L_paper1(i) / C_D2(i)))^2)) * L_D_alpha(i);
%     for j = 1:length(mu)
%         C_M_paper1(j,i) = C_Mc4 - (R_dp/c) * (C_Ds * cos(alpha(i) + mu(j))) ...
%             - R_dp * n * R_dp * d * (cos(alpha(i) + mu(j))^2) / (S * 2 * c) ...
%             - m_s * g * R_dp * sin(alpha(i) + mu(j) - gamma(:,i)) / (0.5 * rho * V^2 * S * c);
% 
%         [~, idx] = min(abs(C_M_paper1(:,i)));
%         corr_mu(i) = mu(idx) * 180/pi;
%         corr_mu(i) = corr_mu(i) * pi/180;
%         C_M2_alpha(i) = (R_dp/c) * (C_Ds * sin(alpha(i) + corr_mu(i))) ...
%             + R_dp * n * R_dp * d * (sin(alpha(i) + corr_mu(i))^2) / (S * 2 * c) ...
%             - d_gamma_alpha * m_s * g * R_dp * cos(alpha(i) + corr_mu(i) - gamma(:,i)) / (0.5 * rho * V^2 * S * c);
%         % C_M2_alpha(i) = C_M2_alpha(i) * pi/180;
%         % Penalize if C_M2_alpha is not negative
%         if C_M2_alpha >= 0
%             corr_mu(i) = Inf;
%         end
%     end
% end
% 
% %updating C_M
% for i =1:length(alpha)
%     up_C_M2(i) = C_Mc4 - (R_dp/c) * (C_Ds * cos(alpha(i) + corr_mu(i))) ...
%         - R_dp * n * R_dp * d * (cos(alpha(i) + corr_mu(i))^2) / (S * 2 * c) ...
%         - m_s * g * R_dp * sin(alpha(i) + corr_mu(i) - gamma(:,i)) / (0.5 * rho * V^2 * S * c);
% end
% Plotting updated C_M2 vs alpha

%% data storing

% headers = {'alpha (deg)', 'C_L2', 'C_L2alpha', 'C_D2', 'C_D2alpha', 'C_L2/C_D2','up_C_M2', 'C_M2_alpha','corr_mu'};
% data = [alpha' .* 180/pi, C_L2(:,2), C_L2alpha(:,2), C_D2(:,2), C_L2(:,2) ./ C_D2(:,2), C_D2alpha(:,2), up_C_M2', C_M2_alpha',corr_mu' .*180/pi];
% filename = 'aerodynamic_analysis_results.xlsx';
% writecell(headers, filename, 'Sheet', 1, 'Range', 'A1');
% writematrix(data, filename, 'Sheet', 1, 'Range', 'A2');
