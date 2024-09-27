% a nicer version of atmosisa
% var = return variable
% height = altitude in m
% idx = atmosisa variable index
%   1 = temperature
%   2 = speed of sound
%   3 = pressure
%   4 = air density

function var = atmos(height, idx) 
	[A, B, C, D] = atmosisa(height);
	mat = [A; B; C; D];
	var = mat(idx, :);
end
