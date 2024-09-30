% Parafoil geometry
% started 26/09/24 - Rosalind Aves

function [b, c, S, AR, t, mu, eps, a, R, d, n, m_s, m_p, A_cube] = calcPfoilGeometry()

b = 1.267; %wing span
c = (0.192 + 0.55)/2; %wing chord

S=b*c; %wing area
t=0.11*c; % airfoil thickness; 
AR = b/c;

mu = deg2rad(10); % rigging angle in radians (imposed incidence by line lengths)

% Surface area and chord length calculations
d = 1.5e-3; % Line diameter (m)
n = 24;
m_p = 0.2;
m_v = 30;
m_s = m_p + m_v;
A_cube = sqrt(0.3^2 + 0.1^2)*0.1; % Cubesat area (m^2)

%% anhedral angle calculations
controlLengthWT = 0.365;

maxRadius = 1.39 - controlLengthWT; % vertical distance from centre of parafoil to end of bridle line

bridle1 = mean([685 702 678 690])/1000;
bridle3 = mean([821 819 837])/1000 ;
bridle5 = mean([870 886 887 905])/1000 ;
bridle7 = mean([955 940 959])/1000 ;

R = mean([bridle1 bridle3 bridle5 bridle7 maxRadius]);

cellLength = b/2 / 3.5; % 3 and a half double cells on each side of parafoil

eps13 = acos((bridle1^2 + bridle3^2 - (cellLength)^2) / (2*bridle1*bridle3));
eps35 = acos((bridle3^2 + bridle5^2 - (cellLength)^2) / (2*bridle3*bridle5));
eps57 = acos((bridle5^2 + bridle7^2 - (cellLength)^2) / (2*bridle5*bridle7));
eps7root = acos((bridle7^2 + maxRadius^2 - (cellLength/2)^2) / (2*bridle7*maxRadius));

eps = sum([eps13 eps35 eps57 eps7root]);

a = b/2 * (1-cos(eps)) / eps; % distance from parafoil tip to root (height) in m
end