function xdot = six_dof_parachute(x, delta, aeroParams, pfoilParams, g)
% format:
% x(1)=X (North-East-Up format)
% x(2)=Y
% x(3)=Z
% x(4)=roll (phi)
% x(5)=pitch (theta)
% x(6)=yaw (psi)
% x(7)=d/dt X
% x(8)=d/dt Y
% x(9)=d/dt Z
% x(10)=d/dt roll
% x(11)=d/dt pitch
% x(12)=d/dt yaw
%
% u(1)=brake
% u(2)=aileron
% xdot=zeros(12,1);
%simpler variables

X = x(1); 
Y = x(2); 
Z = x(3); %X Y Z

phi = x(4);
theta = x(5); 
psi = x(6); %roll pitch yaw

u = x(7); 
v = x(8); 
w = x(9); %U V W

p = x(10); 
q = x(11); 
r = x(12); %P Q R

if (Z < 2000) && (Z > 1000)
    deltaS = delta(1);
    deltaA = delta(2);
else
    deltaS = 0;
    deltaA = 0;
end

m_pay = pfoilParams.m_s - pfoilParams.m_p;
m_p = pfoilParams.m_p;
m = (m_p+m_pay);
mu = pfoilParams.mu;

rho = atmos(Z, 4); %density of air (function of height)
S = pfoilParams.S;
b = pfoilParams.b;
c = pfoilParams.c;
AR = pfoilParams.AR;
eps = pfoilParams.eps;
R = pfoilParams.R;
l_cont = pfoilParams.l_cont;

R_BN = [cos(psi)*cos(theta) sin(psi)*cos(theta) -sin(theta); ...
        cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi) sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi) cos(theta)*sin(phi); ...
        cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi) sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi) cos(theta)*cos(phi)];

R_PB = [cos(mu) 0 -sin(mu); 0 1 0; sin(mu) 0 cos(mu)];

% wind - inertial
W = [0; 0; 0];

% ground speed - body
V = [u; v; w];

% airspeed velocity - body
Va = V - R_BN*W;

% canopy velocity
Vp = R_PB' * V; %R_PB*R_BN'*V

% relative airspeed velocity in canopy frame
Vap = R_PB' * Va;

% local angle of attack
alpha = atan(Vp(3) / Vp(1)); % should be in canopy frame (* R_PB*R_BN'*V ?)

% local sideslip angle
beta = atan(Vp(2) / sqrt(Vp(1)^2 + Vp(3)^2));

% flight path angle
gamma = asin(Va(3) / norm(Va));


R_WB = [cos(alpha)*cos(beta) -sin(beta) sin(alpha)*cos(beta); ...
        cos(alpha)*sin(beta) cos(beta) sin(alpha)*sin(beta); ...
        -sin(alpha) 0 cos(alpha)];

% inertias
Ix_pay = m_pay/12 * (0.1^2 + 0.3^2);
Iy_pay = m_pay/12 * (0.1^2 + 0.3^2);
Iz_pay = m_pay/12 * (0.1^2 + 0.1^2);

v_vol = 0.09 * pfoilParams.c^2 * pfoilParams.b;
% b_inf = 2 * pfoilParams.R * sin(pfoilParams.eps);
h_mean = v_vol / (pfoilParams.c * b);

Ix_p = (v_vol * rho + m_p)/12 * (b^2 + h_mean^2);
Iy_p = (v_vol * rho + m_p)/12 * (c^2 + h_mean^2);
Iz_p = (v_vol * rho + m_p)/12 * (b^2 + c^2);

I = [Ix_p*cos(mu)^2+Iz_p*sin(mu)^2 0 0.5*(Iz_p - Ix_p)*sin(2*mu); ...
     0 Iy_p 0; ...
     0.5*(Iz_p - Ix_p)*sin(2*mu) 0 Ix_p*cos(mu)^2+Iz_p*sin(mu)^2];

% apparent masses and inertias
abar = pfoilParams.a / b;
tbar = pfoilParams.t / c;

A = 0.666 * rho * (1+8/3*abar^2) * pfoilParams.t^2*pfoilParams.b;
B = 0.267 * rho * (1 + 2*abar^2 / tbar^2 * AR^2 * (1-tbar^2)) * pfoilParams.t^2 * c;
C = 0.785 * rho * sqrt(1 + 2*abar^2 * (1-tbar^2)) * AR / (1 + AR) * c^2 * b;

I_A = 0.055 * rho * AR / (1 + AR) * c^2 * b^3;
I_B = 0.0308 * rho * AR / (1 + AR) * (1 + (pi/6) * (1 + AR) * AR * abar^2 * tbar^2) * c^4 * b;
I_C = 0.0555 * rho * (1+8*abar^2) * pfoilParams.t^2 * b^3;

I_am = [A 0 0; 0 B 0; 0 0 C];
% K_abc = eye(3); % AS LESS THAN 1, SO FLOATS
I_ai = [I_A 0 0; 0 I_B 0; 0 0 I_C];
% K_Iabc = eye(3); % AS LESS THAN 1, SO FLOATS

I_am_star = R_PB' * I_am * R_PB;
I_ai_star = R_PB' * I_ai * R_PB;


M_BA = [cos(mu + alpha) 0 -sin(mu + alpha); ...
        0 1 0; ...
        sin(mu + alpha) 0 cos(mu + alpha)];

% Aerodynamic and control coefficients
CL_0= aeroParams.CL0;
CL_alpha= aeroParams.CLalpha;
CL_deltas = aeroParams.CLdeltas;

CD_0= aeroParams.CD0;
CD_s = aeroParams.CDs;
CD_l = (pfoilParams.n .* R .* pfoilParams.d .* cos(alpha).^3) ./ S;
CD_deltas = ((CL_alpha^2/(aeroParams.e*AR*pi)) * (alpha - aeroParams.alpha_zl - aeroParams.dalpha_zl)^2 + aeroParams.dC_D0del) * 2 * aeroParams.bkb;

Cm_q = aeroParams.Cmq;

Cn_p = aeroParams.Cnp;
Cn_r = aeroParams.Cnr; % always negative for stability
Cn_beta = aeroParams.Cnbeta;
% Cn_alphap = aeroParams.Cnalphap; %IGNORING FOR NOW
% Cn_alphar = aeroParams.Cnalphar; %IGNORING FOR NOW
% Cn_alphabeta = aeroParams.Cnalphabeta; %IGNORING FOR NOW
Cn_deltaa = CD_deltas * aeroParams.lk/(2*b);

Cl_p = aeroParams.Clp;
Cl_r = aeroParams.Clr;
Cl_beta = aeroParams.Clbeta;
% Cl_alphar = aeroParams.Clalphar; %IGNORING FOR NOW
Cl_deltaa = aeroParams.Cldeltaa;

Cy_p = aeroParams.Cyp;
Cy_r = aeroParams.Cyr;
Cy_beta = aeroParams.Cybeta;
% Cy_alphar = aeroParams.Cyalphar; %IGNORING FOR NOW
Cy_deltaa = aeroParams.Cydeltaa;


%Coefficient buildup
CL=CL_0 + CL_alpha * alpha + CL_deltas * (2*deltaS + abs(deltaA)); %lift coefficient

CD = CD_0 + CD_s + CD_l + CL_alpha.^2 .* (alpha*cos(eps/2) - aeroParams.alpha_zl).^2 / (aeroParams.e*AR*pi) + CD_deltas*(2*deltaS + abs(deltaA));

CY = Cy_beta * beta + Cy_p * ((b*p)/(2*norm(Vap))) + Cy_r * ((b*r)/(2*norm(Vap))) + Cy_deltaa * deltaA;

Cm = ((c*q)/(2*norm(Vap)))*Cm_q;

Cn = Cn_beta * beta + Cn_p * ((b*p)/(2*norm(Vap))) + Cn_r * ((b*r)/(2*norm(Vap))) + Cn_deltaa * deltaA;

Cl = Cl_beta * beta + Cl_p * ((b*p)/(2*norm(Vap))) + Cl_r * ((b*r)/(2*norm(Vap))) + Cl_deltaa * deltaA;


% Forces
Q = 0.5 * rho * norm(Va)^2;
Fg = m * g * [-sin(theta); cos(theta)*sin(phi); cos(theta)*cos(phi)];
Fa = -Q * S * R_WB' * [CD; CY; CL];

% Moments
Ma = Q * S * [b * Cl; c * Cm; b * Cn];

% Skew symmetric matrices
S_omega = [ 0 -r q; r 0 -p; -q p 0];
Rbm = [0 0 -(l_cont+R)*cos(eps)];
S_rbm = [0 -Rbm(3) 0; ...
         Rbm(3) 0 0; ...
         0 0 0];

I_3x3 = eye(3);


% SOLVE
me = v_vol * rho;

A = [(m + me)*I+I_am_star -I_am_star*S_rbm; ...
      S_rbm*I_am_star (I+I_ai_star-S_rbm'*I_am_star*S_rbm)]; % should A(4) have transpose

% deter = (((m + me)*I+I_am_star)*(I+I_ai_star-S_rbm*I_am_star*S_rbm) - (-I_am_star*S_rbm)*(S_rbm*I_am_star))^-1;


B = [Fa + Fg - S_omega*((m+me)*I_3x3+I_am_star)*[u; v; w] + S_omega*I_am_star*S_rbm*[p; q; r] + S_omega*I_am_star*R_BN*W; ...
     Ma - (S_omega*(I+I_ai_star) - S_rbm*S_omega*I_am_star*S_rbm)*[p; q; r] - S_rbm*S_omega*I_am_star*[u; v; w] + S_rbm*S_omega*I_am_star*R_BN*W];

omegaDOT = A\B; % [udot vdot wdot pdot qdot rdot]'

xdot(1:3) = [1 0 0; 0 1 0; 0 0 -1]* R_BN' * [u; v; w]; % navigation equations (translational)
xdot(4:6) = [1 sin(phi)*tan(theta) cos(phi)*tan(theta); ...
             0 cos(phi) -sin(phi); ...
             0 sin(phi)/cos(theta) cos(phi)/cos(theta)] * [p; q; r]; % kinematic equations (rotational)
xdot(7:9) = omegaDOT(1:3);
xdot(10:12) = omegaDOT(4:6);
% xdot = xdot';


%%
% H=[1 tan(PHI(2))*sin(PHI(1)) tan(PHI(2))*cos(PHI(1)); ...
%     0 cos(PHI(1)) -sin(PHI(1)); ...
%     0 sin(PHI(1))/cos(PHI(2)) cos(PHI(1))/cos(PHI(2))];
% H_inv = H^-1;
% 
% %skew symm cross product matrix
% OMEGA = [ 0 -x(12) x(11); x(12) 0 -x(10); -x(11) x(10) 0];
% Vair = [0 -x(9) x(8); x(9) 0 -x(7); -x(8) x(7) 0];
% 
% %direction cosine matrix
% DCM = [1 0 0; 0 cos(PHI(1)) sin(PHI(1)); ...
%     0 -sin(PHI(1)) cos(PHI(1))] * ...
%     [cos(PHI(2)) 0 -sin(PHI(2)); 0 1 0; ...
%     sin(PHI(2)) 0 cos(PHI(2))] * ...
%     [cos(PHI(3)) sin(PHI(3)) 0; ...
%     -sin(PHI(3)) cos(PHI(3)) 0; 0 0 1];
% 
% % DCM = [cos(PHI(2))*cos(PHI(3)) cos(PHI(2))*sin(PHI(3)) -sin(PHI(2)); ...
% %     (-cos(PHI(1))*sin(PHI(3)) + sin(PHI(1))*sin(PHI(2))*cos(PHI(3))) (cos(PHI(1))*cos(PHI(3)) + sin(PHI(1))*sin(PHI(2))*sin(PHI(3))) sin(PHI(1))*cos(PHI(2)); ...
% %     (sin(PHI(1))*sin(PHI(3)) + cos(PHI(1))*sin(PHI(2))*cos(PHI(3))) (-sin(PHI(1))*cos(PHI(3)) + cos(PHI(1))*sin(PHI(2))*sin(PHI(3))) cos(PHI(1))*cos(PHI(2))];
% 
% 
% v_wind = [0; 0; 0]; % body coordinates
% 
% V = V - v_wind;
% % if sign(V(2)/V(1)) == -1
%     % beta = tan(pi/2 - abs(V(2)/V(1))); % sideslip angle
% % else
%     beta = atan2(V(2), V(1)); % sideslip angle
% % end
% if V(1) == 0 && V(2) == 0
%     beta = 0;
% 
% elseif V(1) == 0 % which direction
%     if V(2) > 0
%         beta = pi/2;
%     elseif V(2) < 0
%         beta = -pi/2;
%     end
% elseif V(2) == 0
%     beta = 0;
% end
% 
% rho = atmos(x(3), 4); %density of air (function of height)
% 
% Xcg = 0; %dist from confluence point to systems cg
% %dist from confluence point to quarter chord
% %point on parafoil along z axis
% Zcg = 1.2;
% 
% rigging_angle = pfoilParams.mu;
% alpha = (atan2(V(3), V(1))) - rigging_angle; % angle of attack 
% 
% if V(1) == 0 && V(3) == 0
%     alpha = 0;
% elseif V(1) == 0 % which direction
%     if V(3) > 0
%         alpha = pi/2;
%     elseif V(3) < 0
%         alpha = -pi/2;
%     end
% elseif V(3) == 0
%     alpha = 0;
% end
% 
% Vr = norm(V); % resultant velocity
% qbar = 0.5 * rho * Vr^2; % dynamic pressure
% mass_t = pfoilParams.m_s; % total system mass in kg
% mass_p = pfoilParams.m_p; % payload mass
% mass_v = mass_t - mass_p; % vehicle mass (3U)
% 
% 
% %apparent mass terms
% % A2 = 0.899 * pi() / 4 * pfoilParams.t^2 * pfoilParams.b; %OLD
% % B2 = 0.39 * pi() / 4 * (pfoilParams.t^2 + 2 * (alpha)^2) * pfoilParams.c; %OLD
% % C2 = 0.783 * pi() / 4 * pfoilParams.c^2 * pfoilParams.b; %OLD
% % 
% % I_A2 = 0.63 * pi() / 48 * pfoilParams.c^2 * pfoilParams.b^3; %OLD
% % I_B2 = 0.817 * 4 / (48*pi()) * pfoilParams.c^4*pfoilParams.b; %OLD
% % I_C2 = 1.001 * pi() / 48 * pfoilParams.t^2 * pfoilParams.b^3; %OLD
% 
% %NEW
% abar = pfoilParams.a / pfoilParams.b;
% tbar = pfoilParams.t / pfoilParams.c;
% 
% A = 0.666 * rho * (1+8/3*abar^2) * pfoilParams.t^2*pfoilParams.b;
% B = 0.267 * rho * (1 + 2*abar^2 / tbar^2 * pfoilParams.AR^2 * (1-tbar^2)) * pfoilParams.t^2 * pfoilParams.c;
% C = 0.785 * rho * sqrt(1 + 2*abar^2 * (1-tbar^2)) * pfoilParams.AR / (1 + pfoilParams.AR) * pfoilParams.c^2 * pfoilParams.b;
% 
% I_A = 0.055 * rho * pfoilParams.AR / (1 + pfoilParams.AR) * pfoilParams.c^2 * pfoilParams.b^3;
% I_B = 0.0308 * rho * pfoilParams.AR / (1 + pfoilParams.AR) * (1 + (pi/6) * (1 + pfoilParams.AR) * pfoilParams.AR * abar^2 * tbar^2) * pfoilParams.c^4 * pfoilParams.b;
% I_C = 0.0555 * rho * (1+8*abar^2) * pfoilParams.t^2 * pfoilParams.b^3;
% 
% K_abc = [A 0 0; 0 B 0; 0 0 C];
% K_abc = eye(3); % AS LESS THAN 1, SO FLOATS
% K_Iabc = [I_A 0 0; 0 I_B 0; 0 0 I_C];
% K_Iabc = eye(3); % AS LESS THAN 1, SO FLOATS
% 
% %inertia matrix
% % J_1=1.3558*[.1357 0 -.0025; 0 .1506 0; -.0025 0 .0203]; %ft lbs force to Joule from masters thesis
% % J_1=1.3558*[6.1298 0 -0.0025; 0 6.150 0; -0.0025 0 .0203]; %ft lbs force to Joule
% % J_1 = [0.42, 0, 0.03; 0, 0.4, 0; 0.03, 0, 0.053]; % from paper
% 
% v_vol = 0.09 * pfoilParams.c^2 * pfoilParams.b;
% b_inf = 2 * pfoilParams.R * sin(pfoilParams.eps);
% h_mean = v_vol / (pfoilParams.c * b_inf);
% 
% Jxxp = (v_vol * rho + 0.2)/12 * (b_inf^2 + h_mean^2);
% Jyyp = (v_vol * rho + 0.2)/12 * (pfoilParams.c^2 + h_mean^2);
% Jzzp = (v_vol * rho + 0.2)/12 * (b_inf^2 + pfoilParams.c^2);
% 
% Jxxv = 3/12 * (0.1^2 + 0.3^2);
% Jyyv = 3/12 * (0.1^2 + 0.3^2);
% Jzzv = 3/12 * (0.1^2 + 0.1^2);
% 
% J_1 = [Jxxv+Jxxp 0 0; 0 Jyyv+Jyyp 0; 0 0 Jzzv+Jzzp];
% 
% % J_1 = [1.622e6 1794.579 1821.736; 1794.579 1.177e7 -68077.582; 1821.736 -68077.582 1.22e7].* 1e-9; %from CAD
% 
% %J_2= J_1+[I A 0 0; 0 I B 0; 0 0 I C];
% %J = J_2 + [0 0 Xcg*Zcg; 0 sqrt(Xcg*Xcg+Zcg*Zcg) 0;
% %Xcg*Zcg 0 0] * (mass+B);
% J = J_1 + mass_t * [Zcg^2 0 0; 0 Zcg^2 0; 0 0 0];
% 
% % Aerodynamic Coefficients - Iacomini & Cerimele
% CL_0= aeroParams.CL0; % 0.091;
% CL_alpha= aeroParams.CLalpha; %0.9;
% CL_deltas = aeroParams.CLdeltas;
% 
% CD_0= aeroParams.CD0; %0.25;
% CD_s = aeroParams.CDs;
% % CD_alpha=0.12;
% CD_l = (pfoilParams.n .* pfoilParams.R .* pfoilParams.d .* cos(alpha).^3) ./ pfoilParams.S;
% CD_deltas = ((CL_alpha^2/(aeroParams.e*pfoilParams.AR*pi)) * (alpha - aeroParams.alpha_zl - aeroParams.dalpha_zl)^2 + aeroParams.dC_D0del) * 2 * aeroParams.bkb;
% 
% % Cm_o=0.35;
% % Cm_alpha= -0.72;
% Cm_q = aeroParams.Cmq;
% 
% Cn_p = aeroParams.Cnp;
% Cn_r = aeroParams.Cnr; %-0.09; %-0.27; % always negative for stability
% Cn_beta = aeroParams.Cnbeta;
% % Cn_alphap = aeroParams.Cnalphap; %IGNORING FOR NOW
% % Cn_alphar = aeroParams.Cnalphar; %IGNORING FOR NOW
% % Cn_alphabeta = aeroParams.Cnalphabeta; %IGNORING FOR NOW
% Cn_deltaa = CD_deltas * aeroParams.lk/(2*pfoilParams.b);
% 
% Cl_p = aeroParams.Clp;
% % Cl_phi = aeroParams.Clphi; % VALUE IS 0
% Cl_r = aeroParams.Clr;
% Cl_beta = aeroParams.Clbeta;
% % Cl_alphar = aeroParams.Clalphar; %IGNORING FOR NOW
% Cl_deltaa = aeroParams.Cldeltaa;
% 
% Cy_p = aeroParams.Cyp;
% Cy_r = aeroParams.Cyr;
% Cy_beta = aeroParams.Cybeta;
% % Cy_alphar = aeroParams.Cyalphar; %IGNORING FOR NOW
% Cy_deltaa = aeroParams.Cydeltaa;
% 
% 
% %Coefficient buildup
% 
% CL=CL_0 + CL_alpha * alpha + CL_deltas * (2*deltaS + abs(deltaA)) %lift coefficient
% 
% % CD=CD_o + CD_alpha * alpha + CD_brake_squared * brake^2 + CD_brake * brake;
% CD = CD_0 + CD_s + CD_l + CL_alpha.^2 .* (alpha*cos(pfoilParams.eps/2) - aeroParams.alpha_zl).^2 / (aeroParams.e*pfoilParams.AR*pi) + CD_deltas*(2*deltaS + abs(deltaA))
% 
% CY = Cy_beta * beta + Cy_p * ((pfoilParams.b*omega(1))/(2*Vr)) + Cy_r * ((pfoilParams.b*omega(3))/(2*Vr)) + Cy_deltaa * deltaA
% 
% Cm = ((pfoilParams.c*omega(2))/(2*Vr))*Cm_q % +Cm_o + Cm_alpha * alpha
% 
% Cn = Cn_beta * beta + Cn_p * ((pfoilParams.b*omega(1))/(2*Vr)) + Cn_r * ((pfoilParams.b*omega(3))/(2*Vr)) + Cn_deltaa * deltaA
% % Cn=Cn_r*b/(2*vc) * (omega(3)) + Cn_aileron*aileron;
% 
% Cl = Cl_beta * beta + Cl_p * ((pfoilParams.b*omega(1))/(2*Vr)) + Cl_r * ((pfoilParams.b*omega(3))/(2*Vr)) + Cl_deltaa * deltaA
% 
% 
% % %% BREAK IF BROKEN
% % if abs(alpha) > deg2rad(15)
% %     warning('Alpha outside 15 degrees');
% % end
% % %%
% 
% %Force & Moment buildup
% F=0.5*rho*pfoilParams.S*Vr^2*(CL*[sin(rigging_angle+alpha); 0; -cos(rigging_angle+alpha)] ...
%     - CD*[cos(rigging_angle+alpha); 0; -sin(rigging_angle+alpha)] + [0; CY; 0]);
% 
% M = [qbar*pfoilParams.S*pfoilParams.b*Cl - mass_t * abs(g)*Zcg*sin(PHI(1)) * cos(PHI(2)); ...
%     qbar*pfoilParams.S*pfoilParams.c*Cm - mass_t * abs(g) * Zcg * sin(PHI(2)); ...
%     qbar*pfoilParams.S*pfoilParams.b*Cn]; %rolling, pitching & yawing moment
% 
% 
% %Navigation Equations
% % xdot(1:3) =  DCM' * V; %[1 0 0; 0 1 0; 0 0 -1]*
% xdot(1:3) =  [1 0 0; 0 1 0; 0 0 -1] * DCM' * x(7:9); %[1 0 0; 0 1 0; 0 0 -1]
% 
% %Force equations
% %xdot(7:9)=[1/(mass+A); 1/(mass+B); 1/(mass+C)] .*
% %(F -cross(omega,V.*[mass+A;mass+B;mass+C])+mass*g*H inv(:,3));
% 
% Faero = 0.5*rho*pfoilParams.S*Vr^2*(CL*[sin(rigging_angle+alpha); 0; -cos(rigging_angle+alpha)] ...
%     - CD*[cos(rigging_angle+alpha); 0; -sin(rigging_angle+alpha)] + [0; CY; 0]);
% Fg = mass_t .* DCM * [0; 0; g];
% Fapp = -K_abc * x(7:9) - cross(x(10:12), K_abc * x(7:9));
% 
% xdot(7:9) = 1/mass_t*(Faero + Fg + Fapp) - OMEGA * x(7:9);
% 
% xdot(7:9)=(eye(3) + K_abc / mass_t)^-1 * ...
%     ((Faero / mass_t) + DCM * [0; 0; g] - cross(omega, x(7:9) + mass_t^-1 * K_abc * x(7:9))); %+g*H_inv(:,3));
% 
% 
% %Kinematic Equations
% xdot(4:6)= H * x(10:12);
% 
% %Moment Equations
% % xdot(10:12)=(J^-1) * (M - OMEGA * J * ([I A; I B; I C]
% %.* omega)+cross(V,V.*[mass+A;mass+B;mass+C]));
% % xdot(10:12)=(J + K Iabc)^-1 * (M - OMEGA * J * omega);
% 
% Mapp = -K_Iabc * x(10:12) - cross(x(4:6), K_Iabc * x(4:6)) - cross(x(7:9), K_abc * x(7:9))
% Mapp2 = -[I_A 0 0; 0 I_B 0; 0 0 I_C] * x(10:12) - OMEGA *  [I_A 0 0; 0 I_B 0; 0 0 I_C] * x(4:6) - Vair* [A 0 0; 0 B 0; 0 0 C] * x(7:9)
% Maero = 0.5 * rho * Vr^2 * pfoilParams.S .* [pfoilParams.b*Cl; pfoilParams.c*Cm; pfoilParams.b*Cn]
% Mg = [-mass_t*g*Zcg*sin(x(5))*cos(x(4)); -mass_t*g*Zcg*sin(x(4)); 0]
% Msum = Mapp + Maero + Mg;
% 
% xdot(10:12) = J^-1 * (Mapp + Maero + Mg - cross(omega, J*omega));
% 
% % xdot(10:12) = (J + K_Iabc)^-1 * ...
% %     (Maero - cross(omega, K_Iabc * omega) - ...
% %     cross(x(7:9),K_abc * x(7:9)) - cross(omega, J*omega)); %OMEGA * J * omega); %
% 
xdot = xdot';
end
