function xdot = six_dof_parachute(x, u, aeroParams, pfoilParams, g)
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

p = [x(1); x(2); x(3)]; %X Y Z
PHI = ([x(4); x(5); x(6)]); %roll pitch yaw
V = [x(7); x(8); x(9)]; %U V W
omega = [x(10); x(11); x(12)]; %P Q R
brake = u(1);
aileron=u(2);

H=[1 tan(PHI(2))*sin(PHI(1)) tan(PHI(2))*cos(PHI(1)); ...
    0 cos(PHI(1)) -sin(PHI(1)); ...
    0 sin(PHI(1))/cos(PHI(2)) cos(PHI(1))/cos(PHI(2))];
H_inv = H^-1;

%skew symm cross product matrix
OMEGA = [ 0 -x(12) x(11); x(12) 0 -x(10); -x(11) x(10) 0];

%direction cosine matrix
DCM = [1 0 0; 0 cos(PHI(1)) sin(PHI(1)); ...
    0 -sin(PHI(1)) cos(PHI(1))] * ...
    [cos(PHI(2)) 0 -sin(PHI(2)); 0 1 0; ...
    sin(PHI(2)) 0 cos(PHI(2))] * ...
    [cos(PHI(3)) sin(PHI(3)) 0; ...
    -sin(PHI(3)) cos(PHI(3)) 0; 0 0 1];

% DCM = [cos(PHI(2))*cos(PHI(3)) cos(PHI(2))*sin(PHI(3)) -sin(PHI(2)); ...
%     (-cos(PHI(1))*sin(PHI(3)) + sin(PHI(1))*sin(PHI(2))*cos(PHI(3))) (cos(PHI(1))*cos(PHI(3)) + sin(PHI(1))*sin(PHI(2))*sin(PHI(3))) sin(PHI(1))*cos(PHI(2)); ...
%     (sin(PHI(1))*sin(PHI(3)) + cos(PHI(1))*sin(PHI(2))*cos(PHI(3))) (-sin(PHI(1))*cos(PHI(3)) + cos(PHI(1))*sin(PHI(2))*sin(PHI(3))) cos(PHI(1))*cos(PHI(2))];


v_wind = [10; 0; 0]; % body coordinates

V = V - v_wind;
% if sign(V(2)/V(1)) == -1
    % beta = tan(pi/2 - abs(V(2)/V(1))); % sideslip angle
% else
    beta = atan2(V(2), V(1)); % sideslip angle
% end
if V(1) == 0 && V(2) == 0
    beta = 0;

elseif V(1) == 0 % which direction
    if V(2) > 0
        beta = pi/2;
    elseif V(2) < 0
        beta = -pi/2;
    end
elseif V(2) == 0
    beta = 0;
end

rho = atmos(x(3), 4); %density of air (function of height)

Xcg = 0; %dist from confluence point to systems cg
%dist from confluence point to quarter chord
%point on parafoil along z axis
Zcg = 1.2;

rigging_angle = pfoilParams.mu;
alpha = (atan2(V(3), V(1))) - rigging_angle; % angle of attack 

Vr = norm(V); % resultant velocity
qbar = 0.5 * rho * Vr^2; % dynamic pressure
mass_t = pfoilParams.m_s; % total system mass in kg
mass_p = pfoilParams.m_p; % payload mass
mass_v = mass_t - mass_p; % vehicle mass (3U)


%apparent mass terms
% A2 = 0.899 * pi() / 4 * pfoilParams.t^2 * pfoilParams.b; %OLD
% B2 = 0.39 * pi() / 4 * (pfoilParams.t^2 + 2 * (alpha)^2) * pfoilParams.c; %OLD
% C2 = 0.783 * pi() / 4 * pfoilParams.c^2 * pfoilParams.b; %OLD
% 
% I_A2 = 0.63 * pi() / 48 * pfoilParams.c^2 * pfoilParams.b^3; %OLD
% I_B2 = 0.817 * 4 / (48*pi()) * pfoilParams.c^4*pfoilParams.b; %OLD
% I_C2 = 1.001 * pi() / 48 * pfoilParams.t^2 * pfoilParams.b^3; %OLD

%NEW
abar = pfoilParams.a / pfoilParams.b;
tbar = pfoilParams.t / pfoilParams.c;

A = 0.666 * rho * (1+8/3*abar^2) * pfoilParams.t^2*pfoilParams.b;
B = 0.267 * rho * (1 + 2*abar^2 / tbar^2 * pfoilParams.AR * (1-tbar^2)) * pfoilParams.t^2 * pfoilParams.c;
C = 0.785 * rho * sqrt(1 + 2*abar^2 * (1-tbar^2)) * pfoilParams.AR / (1 + pfoilParams.AR) * pfoilParams.c^2 * pfoilParams.b;

I_A = 0.055 * rho * pfoilParams.AR / (1 + pfoilParams.AR) * pfoilParams.c^2 * pfoilParams.b^3;
I_B = 0.0308 * rho * pfoilParams.AR / (1 + pfoilParams.AR) * (1 + (pi/6) * (1 + pfoilParams.AR) * pfoilParams.AR * abar^2 * tbar^2) * pfoilParams.c^4 * pfoilParams.b;
I_C = 0.0555 * rho * (1+8*abar^2) * pfoilParams.t^2 * pfoilParams.b^3;

K_abc = [A 0 0; 0 B 0; 0 0 C];
K_abc = eye(3); % AS LESS THAN 1, SO FLOATS
K_Iabc = [I_A 0 0; 0 I_B 0; 0 0 I_C];
K_Iabc = eye(3); % AS LESS THAN 1, SO FLOATS

%inertia matrix
% J_1=1.3558*[.1357 0 -.0025; 0 .1506 0; -.0025 0 .0203]; %ft lbs force to Joule from masters thesis
% J_1=1.3558*[6.1298 0 -0.0025; 0 6.150 0; -0.0025 0 .0203]; %ft lbs force to Joule
% J_1 = [0.42, 0, 0.03; 0, 0.4, 0; 0.03, 0, 0.053]; % from paper

v_vol = 0.09 * pfoilParams.c^2 * pfoilParams.b;
b_inf = 2 * pfoilParams.R * sin(pfoilParams.eps);
h_mean = v_vol / (pfoilParams.c * b_inf);

Jxxp = (v_vol * rho + 0.2)/12 * (b_inf^2 + h_mean^2);
Jyyp = (v_vol * rho + 0.2)/12 * (pfoilParams.c^2 + h_mean^2);
Jzzp = (v_vol * rho + 0.2)/12 * (b_inf^2 + pfoilParams.c^2);

Jxxv = 3/12 * (0.1^2 + 0.3^2);
Jyyv = 3/12 * (0.1^2 + 0.3^2);
Jzzv = 3/12 * (0.1^2 + 0.1^2);

J_1 = [Jxxv+Jxxp 0 0; 0 Jyyv+Jyyp 0; 0 0 Jzzv+Jzzp];

% J_1 = [1.622e6 1794.579 1821.736; 1794.579 1.177e7 -68077.582; 1821.736 -68077.582 1.22e7].* 1e-9; %from CAD

%J_2= J_1+[I A 0 0; 0 I B 0; 0 0 I C];
%J = J_2 + [0 0 Xcg*Zcg; 0 sqrt(Xcg*Xcg+Zcg*Zcg) 0;
%Xcg*Zcg 0 0] * (mass+B);
J = J_1 + mass_t * [Zcg^2 0 0; 0 Zcg^2 0; 0 0 0];

% Aerodynamic Coefficients - Iacomini & Cerimele
CL_0= aeroParams.CL0; % 0.091;
CL_alpha= aeroParams.CLalpha; %0.9;
% CL_brake=0.377;

CD_0= aeroParams.CD0; %0.25;
CD_s = aeroParams.CDs;
% CD_alpha=0.12;
C_Dl = (pfoilParams.n .* pfoilParams.R .* pfoilParams.d .* cos(alpha).^3) ./ pfoilParams.S;
% CD_brake_squared=0.266;
% CD_brake=0.076;

% Cm_o=0.35;
% Cm_alpha= -0.72;
% Cm_brake_squared=0.16;
% Cm_brake=0.056;
Cm_q = aeroParams.Cmq;

Cn_p = aeroParams.Cnp;
Cn_r = aeroParams.Cnr; %-0.09; %-0.27; % always negative for stability
Cn_beta = aeroParams.Cnbeta;
% Cn_alphap = aeroParams.Cnalphap; %IGNORING FOR NOW
% Cn_alphar = aeroParams.Cnalphar; %IGNORING FOR NOW
% Cn_alphabeta = aeroParams.Cnalphabeta; %IGNORING FOR NOW
% Cn_aileron= 0.007; %0.05;

Cl_p = aeroParams.Clp;
% Cl_phi = aeroParams.Clphi; % VALUE IS 0
Cl_r = aeroParams.Clr;
Cl_beta = aeroParams.Clbeta;
% Cl_alphar = aeroParams.Clalphar; %IGNORING FOR NOW

Cy_p = aeroParams.Cyp;
Cy_r = aeroParams.Cyr;
Cy_beta = aeroParams.Cybeta;
% Cy_alphar = aeroParams.Cyalphar; %IGNORING FOR NOW


%Coefficient buildup

CL=CL_0 + CL_alpha * alpha; %lift coefficient

% CD=CD_o + CD_alpha * alpha + CD_brake_squared * brake^2 + CD_brake * brake;
CD = CD_0 + CD_s + C_Dl + CL_alpha.^2 .* (alpha*cos(pfoilParams.eps/2) - aeroParams.alpha_zl).^2 / (aeroParams.e*pfoilParams.AR*pi);

CY = Cy_beta * beta + Cy_p * ((pfoilParams.b*omega(1))/(2*Vr)) + Cy_r * ((pfoilParams.b*omega(3))/(2*Vr));

Cm = ((pfoilParams.c*omega(2))/(2*Vr))*Cm_q; % +Cm_o + Cm_alpha * alpha

Cn = Cn_beta * beta + Cn_p * ((pfoilParams.b*omega(1))/(2*Vr)) + Cn_r * ((pfoilParams.b*omega(3))/(2*Vr));
% Cn=Cn_r*b/(2*vc) * (omega(3)) + Cn_aileron*aileron;

Cl = Cl_beta * beta + Cl_p * ((pfoilParams.b*omega(1))/(2*Vr)) + Cl_r * ((pfoilParams.b*omega(3))/(2*Vr));

%Force & Moment buildup
F=0.5*rho*pfoilParams.S*Vr^2*(CL*[sin(rigging_angle+alpha); 0; -cos(rigging_angle+alpha)] ...
    - CD*[cos(rigging_angle+alpha); 0; -sin(rigging_angle+alpha)] + [0; CY; 0]);

M = [qbar*pfoilParams.S*pfoilParams.b*Cl - mass_t * abs(g)*Zcg*sin(PHI(1)) * cos(PHI(2)); ...
    qbar*pfoilParams.S*pfoilParams.c*Cm - mass_t * abs(g) * Zcg * sin(PHI(2)); ...
    qbar*pfoilParams.S*pfoilParams.b*Cn]; %rolling, pitching & yawing moment


%Navigation Equations
% xdot(1:3) =  DCM' * V; %[1 0 0; 0 1 0; 0 0 -1]*
xdot(1:3) =  [1 0 0; 0 1 0; 0 0 -1] * DCM' * x(7:9); %[1 0 0; 0 1 0; 0 0 -1]

%Force equations
%xdot(7:9)=[1/(mass+A); 1/(mass+B); 1/(mass+C)] .*
%(F -cross(omega,V.*[mass+A;mass+B;mass+C])+mass*g*H inv(:,3));

Faero = 0.5*rho*pfoilParams.S*Vr^2*(CL*[sin(rigging_angle+alpha); 0; -cos(rigging_angle+alpha)] ...
    - CD*[cos(rigging_angle+alpha); 0; -sin(rigging_angle+alpha)] + [0; CY; 0]);
Fg = mass_t .* DCM * [0; 0; g];
Fapp = -K_abc * x(7:9) - cross(x(10:12), K_abc * x(7:9));

xdot(7:9) = 1/mass_t*(Faero + Fg + Fapp) - OMEGA * x(7:9);
% xdot(7:9)=(eye(3) + K_abc / mass)^-1 * ...
%     ((F / mass) - cross(omega, x(7:9) + mass^-1 * K_abc * x(7:9) + G)); %+g*H_inv(:,3));
% [(eye(3) + K abc / mass); alpha V(3) V(1)]

%Kinematic Equations
xdot(4:6)= H * x(10:12);

%Moment Equations
% xdot(10:12)=(J^-1) * (M - OMEGA * J * ([I A; I B; I C]
%.* omega)+cross(V,V.*[mass+A;mass+B;mass+C]));
% xdot(10:12)=(J + K Iabc)^-1 * (M - OMEGA * J * omega);

Mapp = -K_Iabc * x(10:12) - cross(x(4:6), K_Iabc * x(4:6)) - cross(x(7:9), K_abc * x(7:9));
Maero = 0.5 * rho * Vr^2 * pfoilParams.S .* [pfoilParams.b*Cl; pfoilParams.c*Cm; pfoilParams.b*Cn];
Mg = [-mass_t*g*Zcg*sin(x(5))*cos(x(4)); -mass_t*g*Zcg*sin(x(4)); 0];
Msum = Mapp + Maero + Mg;

xdot(10:12) = J^-1 * (Mapp + Maero + Mg - cross(omega, J*omega));
% xdot(10:12) = (J + K_Iabc)^-1 * ...
%     (M - cross(omega, K_Iabc * omega) - ...
%     cross(V,K_abc * V) - cross(omega, J*omega)) %OMEGA * J * omega); %

xdot = xdot';
end