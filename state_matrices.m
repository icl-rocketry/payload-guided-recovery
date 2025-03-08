syms x1 x2 x3 x4 x5 x6 
syms deltaA % roll angle
syms W0

cell = num2cell(pfoilList);
[pfoilParams.b, pfoilParams.c, pfoilParams.S, pfoilParams.AR, pfoilParams.t, pfoilParams.mu, pfoilParams.eps, pfoilParams.a, pfoilParams.R, pfoilParams.d, pfoilParams.n, pfoilParams.m_s, pfoilParams.m_p, pfoilParams.A_cube, pfoilParams.l_cont] = cell{:};
Va_norm = x1;
gamma = x2;
Psi = x3; % heading angle
% X = x4;
% Y = x5;
% Z = x6;

% if (abs(Z) < 2825) % pretied
%     deltaA = -delta(2);
% else
%     deltaA = 0;
% end

% if (abs(Z) < 2491) && (abs(Z) > 2357) % 50 sec actuation
%     deltaA = delta(2);
% elseif (abs(Z) < 2357) && (abs(Z) > 2241) % 20 sec hold neutral
%     deltaA = 0;
% elseif (abs(Z) < 2241) % rest of flight actuation
%     deltaA = -delta(2);
% else
%     deltaA = 0;
% end

S = pfoilParams.S;
m = pfoilParams.m_s;
rho = 1.225; % atmos(Z, 4);
phi = deltaA; % bank angle (set to 0 if no control)
W = W0; % wind vector - inertial

R_WN = [cos(Psi)*cos(gamma) sin(Psi)*cos(gamma) -sin(gamma); ...
        cos(Psi)*sin(gamma)*sin(phi) - sin(Psi)*cos(phi) sin(Psi)*sin(gamma)*sin(phi) + cos(Psi)*cos(phi) cos(gamma)*sin(phi); ...
        cos(Psi)*sin(gamma)*cos(phi) + sin(Psi)*sin(phi) sin(Psi)*sin(gamma)*cos(phi) - cos(Psi)*sin(phi) cos(gamma)*cos(phi)];

V = R_WN' * [Va_norm; 0; 0] + R_WN * W;
Va = V - R_WN * W;
alpha = atan(abs(Va(3) / Va(1)));


% Aerodynamic and control coefficients
CL_0 = aeroList(3);
CL_alpha= 0.9; %aeroParams.CLalpha;
CL_deltas = aeroList(28);

CD_0= aeroList(5);
CD_s = aeroList(6);
CD_l = (pfoilParams.n .* pfoilParams.R .* pfoilParams.d .* cos(alpha).^3) ./ S;
CD_deltas = ((CL_alpha^2/(aeroList(2)*pfoilParams.AR*pi)) * (alpha - aeroList(1) - aeroList(25))^2 + aeroList(27)) * 2 * aeroList(23);


CL = CL_0 + CL_alpha * alpha + CL_deltas * (abs(deltaA)); %lift coefficient
CD = CD_0 + CD_s + CD_l + CL_alpha.^2 .* (alpha*cos(pfoilParams.eps/2) - aeroList(1)).^2 / (aeroList(2)*pfoilParams.AR*pi) + CD_deltas*(abs(deltaA));

% Lift and Drag
% CL = 0.2985;
% CD = CL;
L = 0.5 * rho * Va_norm^2 * S * CL;
D = 0.5 * rho * Va_norm^2 * S * CD;

% A = [sin(gamma)*Va -cos(gamma)*cos(phi)*Va; (m*cos(gamma)*Va) (m*cos(gamma)*sin(phi)*Va)];
% b = [-g*cos(gamma)*sin(gamma); (L-m*g*cos(gamma)*cos(phi))];

% omegaDOT = A\b;
xdot1 = -1/m * (D + m*g*sin(gamma));
xdot2 = 1/(m*Va_norm) * (L*cos(phi) - m*g*cos(gamma));
xdot3 = (L*sin(phi)) / (m*Va_norm*cos(gamma));
xdot_3x1 = R_WN' * [Va_norm; 0; 0] + W;
xdot4 = xdot_3x1(1);
xdot5 = xdot_3x1(2);
xdot6 = xdot_3x1(3);

A = jacobian([xdot1, xdot2, xdot3, xdot4, xdot5, xdot6], ...
    [x1, x2, x3, x4, x5, x6]);
B = jacobian([xdot1, xdot2, xdot3, xdot4, xdot5, xdot6], ...
    deltaA);

A=double(subs(A, {x1, x2, x3, x4, x5, x6, deltaA}, {10, 0.1, 0, 0, 0, 0, 0.1}));
B=double(subs(B, {x1, x2, x3, x4, x5, x6, deltaA}, {10, 0.1, 0, 0, 0, 0, 0.1}));

disp(null(ctrb(A, B)))

% lqr
Q = eye(6);
R = 1;
[K, S, P] = lqr(A, B, Q, R);
