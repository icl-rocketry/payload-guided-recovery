% 3 DOF from JPL DARTS 
% started 06/10/24 - Rosalind Aves

function xdot = three_dof_parachute(x, delta, W0, aeroParams, pfoilParams, g)

Va_norm = x(1);
gamma = x(2);
Psi = x(3); % heading angle
X = x(4);
Y = x(5);
Z = x(6);

if (abs(Z) < 2568) && (abs(Z) > 2361)
    deltaA = delta(2);
elseif (abs(Z) < 2361) && (abs(Z) > 2200)
    deltaA = 0;
elseif (abs(Z) < 2275)
    deltaA = delta(2);
else
    deltaA = 0;
end

S = pfoilParams.S;
m = pfoilParams.m_s;
rho = atmos(Z, 4);
phi = deltaA; % bank angle (set to 0 if no control)
W = W0; % wind vector - inertial

R_WN = [cos(Psi)*cos(gamma) sin(Psi)*cos(gamma) -sin(gamma); ...
        cos(Psi)*sin(gamma)*sin(phi) - sin(Psi)*cos(phi) sin(Psi)*sin(gamma)*sin(phi) + cos(Psi)*cos(phi) cos(gamma)*sin(phi); ...
        cos(Psi)*sin(gamma)*cos(phi) + sin(Psi)*sin(phi) sin(Psi)*sin(gamma)*cos(phi) - cos(Psi)*sin(phi) cos(gamma)*cos(phi)];

V = R_WN' * [Va_norm; 0; 0] + R_WN * W;
Va = V - R_WN * W;
alpha = atan(abs(Va(3) / Va(1)));


% Aerodynamic and control coefficients
CL_0 = aeroParams.CL0;
CL_alpha= 0.9; %aeroParams.CLalpha;
CL_deltas = aeroParams.CLdeltas;

CD_0= aeroParams.CD0;
CD_s = aeroParams.CDs;
CD_l = (pfoilParams.n .* pfoilParams.R .* pfoilParams.d .* cos(alpha).^3) ./ S;
CD_deltas = ((CL_alpha^2/(aeroParams.e*pfoilParams.AR*pi)) * (alpha - aeroParams.alpha_zl - aeroParams.dalpha_zl)^2 + aeroParams.dC_D0del) * 2 * aeroParams.bkb;


CL = CL_0 + CL_alpha * alpha + CL_deltas * (abs(deltaA)); %lift coefficient
CD = CD_0 + CD_s + CD_l + CL_alpha.^2 .* (alpha*cos(pfoilParams.eps/2) - aeroParams.alpha_zl).^2 / (aeroParams.e*pfoilParams.AR*pi) + CD_deltas*(abs(deltaA));

% Lift and Drag
% CL = 0.2985;
% CD = CL;
L = 0.5 * rho * Va_norm^2 * S * CL;
D = 0.5 * rho * Va_norm^2 * S * CD;

% A = [sin(gamma)*Va -cos(gamma)*cos(phi)*Va; (m*cos(gamma)*Va) (m*cos(gamma)*sin(phi)*Va)];
% b = [-g*cos(gamma)*sin(gamma); (L-m*g*cos(gamma)*cos(phi))];

% omegaDOT = A\b;

xdot(1,1) = -1/m * (D + m*g*sin(gamma));
xdot(2,1) = 1/(m*Va_norm) * (L*cos(phi) - m*g*cos(gamma));
xdot(3,1) = (L*sin(phi)) / (m*Va_norm*cos(gamma));

xdot(4:6,1) = R_WN' * [Va_norm; 0; 0] + W;
end
