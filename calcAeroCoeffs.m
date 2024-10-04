% CLARK Y Parafoil Aerodynamic Analysis DESIGN POINTS
% started 25/09/24 - Rosalind Aves

function [aeroParams, pfoilParams] = calcAeroCoeffs(pfoilParams, u)
% Constants and parameters
a0 = 0.06*180/pi; % Clark Y lift curve slope (rad^-1)
alpha_zl = -3 * pi/180; % Zero lift angle of attack (rad)

% Span and line length ratio calculations
b = pfoilParams.b; % Span (m)
eps = pfoilParams.eps;
% R = pfoilParams.R;
AR = pfoilParams.AR;
% c = pfoilParams.c;
S = pfoilParams.S;
A_cube = pfoilParams.A_cube;
% mu = pfoilParams.mu;

% beta = b ./ (4 .* R);

%% C_L coefficients !!
% tau = 0.0; % for elliptical planforms
% a0_dash = (2 * pi .* AR) .* tanh(a0 ./ (2 * pi .* AR));
% C_Lalpha_paper1 = (pi .* AR .* a0_dash) ./ (pi .* AR + a0_dash .* (1 + tau)); % FOR AIRFOIL
% C_L2alpha_paper1 = C_Lalpha_paper1.* cos(beta).^2; % FOR PARAFOIL per rad
% 
% C_L_paper1(:) = C_Lalpha_paper1 .* (alpha - alpha_zl) .* cos(beta).^2; %FOR PARAFOIL per rad
% C_L0_paper1 = C_Lalpha_paper1 .* (0 - alpha_zl) .* cos(beta).^2;

k_paper2 = AR*pi/a0;
k1_paper2 = (sqrt(k_paper2^2 + 1) + 1)/(sqrt(k_paper2^2 + 4) + 2);
k2_paper2 = (sqrt(k_paper2^2 + 4) - 1)/(sqrt(k_paper2^2 + 1) + 1);

C_Lalpha_paper2 = AR*pi/(sqrt(k_paper2^2 + 1)+1);
% C_L_paper2(:) = C_Lalpha_paper2.*(alpha.*cos(eps/2)-alpha_zl).*cos(eps/2);
C_L0_paper2 = C_Lalpha_paper2.*(0 .*cos(eps/2)-alpha_zl).*cos(eps/2);


%% C_D coefficients !! - no derivatives just pure CD

% delta = [0.023]; 
e = 0.7;
C_D0_paper1 = 0.015 + 0.004 + 0.5 * 0.117 + 0.0001;
% C_D_paper1 = C_D0_paper1 + ((C_Lalpha_paper1 * (alpha - alpha_zl)).^2 .* (1 + delta)) ./ (pi * AR);

% C_Dl = (n .* R .* d .* cos(alpha).^3) ./ S;
C_Ds = 1.05*A_cube/S; % cuboid CD = 1.05

% C_D_paper1 = C_D_paper1 + C_Dl + C_Ds; % ADDING LINE AND STORE DRAG
% C_D_paper2 = C_D0_paper1 + C_Ds + C_Dl + C_Lalpha_paper2.^2 .* (alpha*cos(eps/2) - alpha_zl).^2 / (e*AR*pi); %MOVE TO 6DOF

%% C_M coefficients !!

% maintaining stable equilibrium at each alpha
% C_Mc4 = -0.08; % C_M0

C_Mq_paper2 = - C_Lalpha_paper2 .* (cos(eps/2)^2)/12;%pitch damping derivative assuming that the CP and reference point coincide

% gamma = atan(1./(C_L_paper1 ./ C_D_paper1));
% 
% Wt = m_s * g;
% Ct = sqrt(C_L_paper1.^2 + C_D_paper1.^2);
% V = sqrt((2 * Wt) ./ (rho .* S .* Ct));
% u = V .* cos(gamma);
% w = V .* sin(gamma);

% C_Ll = - (n * d * R * cos(alpha + mu).^3 .* sin(alpha + mu)) ./ (S);
% C_Ls = 0;
% 
% C_M_paper1 = C_Mc4 - (R/c) * (C_Ds * cos(alpha + mu)) ...
%             - R./ (S * 2 * c) * (n * R * d .* (cos(alpha + mu).^2))  ...
%             - (m_s * g * R * sin(alpha + mu - gamma)) ./ (0.5 .* rho * V.^2 .* S .* c);

%% Yaw damping !!

C_Yr_paper2 = C_Lalpha_paper2 * (sin(eps)/2) * alpha_zl;
C_Yralpha_paper2 = -C_Lalpha_paper2 * (sin(eps)/2) * cos(eps/2)^2;
C_lr_paper2 = -C_Lalpha_paper2 * (sin(eps)/(4*eps)) * alpha_zl;
C_lalphar_paper2 = C_Lalpha_paper2 * (sin(eps/2)/(2*eps)) * cos(eps/2)^2;
C_nr_paper2 = -(C_D0_paper1 * sin(eps))/(3*eps) - (C_Lalpha_paper2^2 * sin(eps/2)^2 * alpha_zl)/(e*AR*pi*8*eps^2) - (C_Lalpha_paper2 * eps^2) / (AR^2 * 24);
C_nalphar_paper2 = ( ((C_Lalpha_paper2^2) / (e * AR * pi)) * (sin(eps)/eps) - C_Lalpha_paper2 * ((sin(eps)^2) / (eps^2)) ) * alpha_zl;

%% Roll damping !!

C_lp_paper2 = - C_Lalpha_paper2 * k1_paper2 * sin(eps)/(8*eps);
C_lphi_paper2 = 0; % no value 

C_Yp_paper2 = C_Lalpha_paper2 * k1_paper2 * sin(eps) / 4;
C_np_paper2 = C_Lalpha_paper2 * k1_paper2 * k2_paper2 * (sin(eps)/(8*eps)) * alpha_zl;
C_nalphap_paper2 = -C_Lalpha_paper2 * k1_paper2 * k2_paper2 * (sin(3*eps/2) / (4*eps));

%% Side slip angle derivatives !!

C_Ybeta_paper2 = -C_Lalpha_paper2 * k1_paper2 * ((eps * sin(eps)) / 4) - C_D0_paper1 * ((1 + 2*cos(eps)) / 3);
C_lbeta_paper2 = C_Lalpha_paper2 * k1_paper2 * (sin(eps) / 8);
C_nbeta_paper2 = -C_Lalpha_paper2 * k1_paper2  * k2_paper2 * (sin(eps) / 8) * alpha_zl;
C_nalphabeta_paper2 = C_Lalpha_paper2 * k1_paper2 * k2_paper2 * ((sin(3*eps/2))/4);

%% Control derivatives (symmetric deflection control has been assumed)
deltaS = u(1); deltaA = u(2); % control inputs
deltaL = (2 * deltaS - deltaA) / 2; % left deflection
deltaR = 2*deltaS - deltaL; % right deflection
b_r = b/2 * cos(deltaR);
b_l = b/2 * cos(deltaL);
b_k = b_r + b_l; % total deflected span 

bk_b = b_k / b; %ratio of deflected surface width to the span assumed to 0.25 for now
 
eps_k_paper2 = eps * (1 - bk_b);
del_alpha_zl_paper2 = -11 * pi/180 ; % (CONSTANT)
lk_paper2 = b/2 * (1 - bk_b)*(sin(eps_k_paper2)/eps_k_paper2);
if eps_k_paper2 == 0
    lk_paper2 = 0;
end
del_C_D0del = 0.2; % additional profile drag taken from the paper (CONSTANT)

C_Ldelta_s_paper2 = -C_Lalpha_paper2 * del_alpha_zl_paper2 * 2 * bk_b * cos(eps_k_paper2);
% C_Ddelta_s_paper2 = (C_Lalpha_paper2^2 * (alpha - alpha_zl - del_alpha_zl)^2/(e*AR*pi) + del_C_D0del) * 2 * bk_b;
C_Ydelta_a_paper2 = C_Ldelta_s_paper2 * sin(eps_k_paper2)/(2*cos(eps_k_paper2));
C_ldelta_a_paper2 = - C_Ldelta_s_paper2/2 * cos(eps_k_paper2/2)/cos(eps_k_paper2)*lk_paper2/b;
% C_ndelta_a_paper2 = C_Ddelta_s_paper2 * lk/(2*b);

%% put into struct
aeroParams.alpha_zl = alpha_zl;
aeroParams.e = e;

aeroParams.CL0 = C_L0_paper2;
aeroParams.CLalpha = C_Lalpha_paper2;

aeroParams.CD0 = C_D0_paper1;
aeroParams.CDs = C_Ds;

aeroParams.Cmq = C_Mq_paper2;

aeroParams.Cnp = C_np_paper2;
aeroParams.Cnr = C_nr_paper2;
aeroParams.Cnbeta = C_nbeta_paper2;
aeroParams.Cnalphap = C_nalphap_paper2;
aeroParams.Cnalphar = C_nalphar_paper2;
aeroParams.Cnalphabeta = C_nalphabeta_paper2;

aeroParams.Clp = C_lp_paper2;
aeroParams.Clphi = C_lphi_paper2;
aeroParams.Clr = C_lr_paper2;
aeroParams.Clbeta = C_lbeta_paper2;
aeroParams.Clalphar = C_lalphar_paper2;

aeroParams.Cyp = C_Yp_paper2;
aeroParams.Cyr = C_Yr_paper2;
aeroParams.Cybeta = C_Ybeta_paper2;
aeroParams.Cyalphar = C_Yralpha_paper2;

aeroParams.bkb = bk_b;
aeroParams.epsk = eps_k_paper2;
aeroParams.dalpha_zl = del_alpha_zl_paper2;
aeroParams.lk = lk_paper2;
aeroParams.dC_D0del = del_C_D0del;
aeroParams.CLdeltas = C_Ldelta_s_paper2;
aeroParams.Cydeltaa = C_Ydelta_a_paper2;
aeroParams.Cldeltaa = C_ldelta_a_paper2;
end