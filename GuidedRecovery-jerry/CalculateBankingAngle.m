%Given chord length l, and distance between chord on parafoil d, calculate
%banking angle phi, given the retraction of chord delta L(dl)
%in cm
l = 90;
d = 100;
theta = acos(0.5*d/l);

dl = 0;

phiArray = [];
dlArray = 0:0.2:20;

for dl = dlArray
    phi = findPhi(l,d,dl);
    disp("dl = " + dl )
    disp("Phi = " + phi*180/pi)
    disp("     ")
    phiArray = [phiArray, real(phi)*180/pi];
end

figure(1)
clf;
plot(dlArray,phiArray,'LineWidth',3)

xlabel("dl [cm]")
ylabel("phi [degrees]")
% Increase the font size of the x-axis label
xlabel_handle = xlabel('dl [cm]');
set(xlabel_handle, 'FontSize', 30);

% Increase the font size of the y-axis label
ylabel_handle = ylabel('phi [degrees]');
set(ylabel_handle, 'FontSize', 30);

ax = gca; % Get current axes
set(ax, 'FontSize', 30); % Set the font size for the tick labels
grid on
% 
% syms phi
% syms c
% 
% eqn1 = c^2 - 2*(l-dl)*cos(phi + theta)*c + (l-dl)^2 - l^2 == 0;
% sol = solve(eqn1,c);
% c = sol(2);
% 
% eqn2 = c^2 + d^2 - 2*c*d*cos(phi) - dl^2 == 0;
% sol = solve(eqn2, phi);
% double(sol)

coefficients = polyfit(dlArray(1:51),phiArray(1:51),1);
xspace = 0:0.1:20;
yspace = coefficients(1)*xspace + coefficients(2);
hold on
plot(xspace,yspace,'LineStyle','--','LineWidth',3)
legend("Actual","phi = 0.5031*dl - 0.0437")
