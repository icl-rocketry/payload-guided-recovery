% Parafoil dimensions
% started 25/09/24 - Rosalind Aves
function [eps, a, R] = calcEpsilon(b, controlLengthWT)

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

% a_forced =  0.25*5*0.33;
% 
% % Newton Raphson
% 
% fun = @(x) (b/2 * (1-cos(x)) / x) - a_forced;
% x = fzero(fun, deg2rad(90));

end