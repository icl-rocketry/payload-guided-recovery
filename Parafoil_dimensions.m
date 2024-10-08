% Parafoil dimensions
% started 25/09/24 - Rosalind Aves

b = 1.27; % m span

bridle1 = mean([685 702 678 690])/1000;
bridle3 = mean([821 819 837])/1000;
bridle5 = mean([870 886 887 905])/1000;
bridle7 = mean([955 940 959])/1000;

cellLength = b/2 / 3.5;
maxRadius = 1.39;
controlLengthWT = 0.365; % line length for maxRadius (function of)

eps13 = acos((bridle1^2 + bridle3^2 - (cellLength)^2) / (2*bridle1*bridle3));
eps35 = acos((bridle3^2 + bridle5^2 - (cellLength)^2) / (2*bridle3*bridle5));
eps57 = acos((bridle5^2 + bridle7^2 - (cellLength)^2) / (2*bridle5*bridle7));


eps7root = acos((bridle5^2 + bridle7^2 - (cellLength/2)^2) / (2*bridle5*bridle7));
eps = sum([eps13 eps35 eps57]);
eps_deg = rad2deg(eps);

a = b/2 * (1-cos(eps)) / eps;