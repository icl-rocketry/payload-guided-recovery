% Wind tunnel results
% started 05/10/24 - Rosalind Aves

S = 0.47;
mu = deg2rad(-5);

F_D = [7,5.5,4]; % drag force from load cell at 9m/s
V = [9,8,7]; % m/s in wind tunnel
rho = 1.225;
alpha = [deg2rad(60), deg2rad(45), deg2rad(30)];

C_D = F_D ./ (0.5 * rho .* V.^2 .* S);

C_N = C_D ./ cos(alpha);
C_L = C_N .* sin(alpha);
F_L = 0.5 .* rho .* S .* C_L .* V.^2;

gamma = atan(C_D ./ C_L);

W = 2.4 * 9.81;

V_wt = sqrt(W./(0.5 .* C_N .* S .* rho));
w_wt = V_wt .* sin(alpha + mu);
u_wt = V_wt .* cos(alpha + mu);

% Mach
w_mach_1 = 6.5; % terminal velocity
gamma_mach = deg2rad(45); % think so ?
V_mach = w_mach_1 / sin(gamma_mach);
W_mach = 9.81;

C_N_mach = W_mach / (0.5 * rho * V_mach^2 * S);

V_mach_3 = sqrt(W/(0.5 * C_N_mach * S * rho));
w_mach_3 = V_mach_3 * sin(gamma_mach);
u_mach_3 = V_mach_3 * cos(gamma_mach);




%% drift and wind

Uwind = 8.7; % max u velocity

time_mach_3 = 3000/w_mach_3;
drift_mach_3 = u_mach_3 * time_mach_3;
drift_mach_3_wind = (u_mach_3 + Uwind) * time_mach_3;

time_wt = 3000./w_wt;
drift_wt = u_wt .* time_wt;
drift_wt_wind = (u_wt + Uwind) .* time_wt;


%% plot CL alpha
CL_a(1) = (C_L(1) - C_L(2)) / (alpha(1) - alpha(2));
CL_a(2) = (C_L(1) - C_L(3)) / (alpha(1) - alpha(3));
CL_a(3) = (C_L(3) - C_L(2)) / (alpha(3) - alpha(2));

CLalpha_wt = mean(CL_a);

% y = mx + CL0
CL0 = C_L(1) - CLalpha_wt*alpha(1);
alpha0 = -CL0/CLalpha_wt;
%  

x = linspace(alpha0, alpha(1));
y = CLalpha_wt.*x + CL0;
figure();
plot(alpha, C_L, 'o')
hold on;
plot(x, y)
