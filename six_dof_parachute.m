function [xdot, y] = six_dof_parachute(x, u)
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
xdot=zeros(12,1);
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

%cross product matrix
OMEGA = [ 0 -x(12) x(11); x(12) 0 -x(10); -x(11) x(10) 0];

%direction cosine matrix
DCM = [1 0 0; 0 cos(PHI(1)) sin(PHI(1)); ...
    0 -sin(PHI(1)) cos(PHI(1))] * ...
    [cos(PHI(2)) 0 sin(PHI(2)); 0 1 0; ...
    -sin(PHI(2)) 0 cos(PHI(2))] * ...
    [cos(PHI(3)) sin(PHI(3)) 0; ...
    -sin(PHI(3)) cos(PHI(3)) 0; 0 0 1];

v_wind = [0; 0; 0];

V = V - v_wind;

%constants
g=9.81; %gravity
rho=1.225; %density of air

a=36 * 0.0254; %height of parafoil
b=1.26; %84 * 0.0254; %wing span 1.26;
c=0.25; %24 * 0.0254; %wing chord 0.25; 
S=b*c; %wing area
t=0.14*c; %5*0.0254; %thickness 0.14*c; 

Xcg = 0; %dist from confluence point to systems cg
%dist from confluence point to quarter chord
%point on parafoil along z axis
Zcg = 48 * 0.0254;

rigging_angle = 10; %degrees

vc = sqrt(V'*V); %velocity
qbar = 0.5 * rho * vc*vc; %dynamic pressure
mass = 3; %4; %mass in kg
alpha = rad2deg(atan2(V(3), V(1))) - rigging_angle;

%apparent mass terms
A = 0.899 * pi() / 4 * t^2 * b;
B = 0.39 * pi() / 4 * (t^2 + 2 * deg2rad(alpha)^2) * c;
C = 0.783 * pi() / 4 * c^2 * b;

I_A = 0.63 * pi() / 48 * c^2 * b^3;
I_B = 0.817 * 4 / (48*pi()) * c^4*b;
I_C = 1.001 * pi() / 48 * t^2 * b^3;


K_abc = [A 0 0; 0 B 0; 0 0 C];
K_Iabc = [I_A 0 0; 0 I_B 0; 0 0 I_C];

[K_abc; K_Iabc];

% K abc = [0 0 0; 0 0 0; 0 0 0];
% K Iabc = [0 0 0; 0 0 0; 0 0 0];

%inertia matrix
J_1=1.3558*[.1357 0 .0025; 0 .1506 0; .0025 0 .0203];

%J_2= J_1+[I A 0 0; 0 I B 0; 0 0 I C];
%J = J_2 + [0 0 Xcg*Zcg; 0 sqrt(Xcg*Xcg+Zcg*Zcg) 0;
%Xcg*Zcg 0 0] * (mass+B);
J = J_1 + mass * [Zcg^2 0 0; 0 Zcg^2 0; 0 0 0];

%Aerodynamic Coefficients - Iacomini & Cerimele
CL_o=0.2;
CL_alpha=0.0375;
CL_brake=0.377;

CD_o=0.12;
CD_alpha=0.0073;
CD_brake_squared=0.266;
CD_brake=0.076;

Cm_o=-0.0115;
Cm_alpha=-0.004;
Cm_brake_squared=0.16;
Cm_brake=0.056;

Cn_r=-0.0936; %always negative for stability
Cn_aileron=0.05;

%Coefficient buildup

CL=CL_o + CL_alpha * alpha + CL_brake*brake; %lift coefficient

CD=CD_o + CD_alpha * alpha + CD_brake_squared * brake^2 + CD_brake * brake;

Cm = Cm_o + Cm_alpha * alpha + Cm_brake*brake + Cm_brake_squared*brake^2;

Cn=Cn_r*b/(2*vc) * rad2deg(omega(3)) + Cn_aileron*aileron;
Cl = 0;

%Force & Moment buildup
F=0.5*rho*S*vc^2*(CL*[sind(rigging_angle+alpha); 0; -cosd(rigging_angle+alpha)] ...
    - CD*[cosd(rigging_angle+alpha); 0; -sind(rigging_angle+alpha)]);

M = [qbar*S*b*Cl - mass * abs(g)*Zcg*sin(PHI(1)) * cos(PHI(2)); ...
    qbar*S*c*Cm - mass * abs(g) * Zcg * sin(PHI(2)); ...
    qbar*S*b*Cn]; %rolling, pitching & yawing moment

[alpha; V;F; vc^2];

%Force equations
%xdot(7:9)=[1/(mass+A); 1/(mass+B); 1/(mass+C)] .*
%(F -cross(omega,V.*[mass+A;mass+B;mass+C])+mass*g*H inv(:,3));
xdot(7:9)=(eye(3) + K_abc / mass)^-1 * ...
    ((F / mass) -cross(omega,V + mass^-1 * K_abc * V)+g*H_inv(:,3));
% [(eye(3) + K abc / mass); alpha V(3) V(1)]

%Kinematic Equations
xdot(4:6)= H * omega;

%Moment Equations
% xdot(10:12)=(J^-1) * (M - OMEGA * J * ([I A; I B; I C]
%.* omega)+cross(V,V.*[mass+A;mass+B;mass+C]));
% xdot(10:12)=(J + K Iabc)^-1 * (M - OMEGA * J * omega);
xdot(10:12)=(J + K_Iabc)^-1 * ...
    (M -cross(omega,K_Iabc * omega) - ...
    cross(V,K_abc * V) - OMEGA * J * omega);


%Navigation Equations
xdot(1:3) = [1 0 0; 0 1 0; 0 0 -1] * DCM' * V;

end