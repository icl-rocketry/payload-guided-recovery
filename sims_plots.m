% trajectory performance
% started 10/10/24 - Rosalind Aves

load drift_control.mat
load drift_nocontrol.mat
load drift_pretie.mat

figure();
p3 = plot3(X_control, Y_control, Z_control);
hold on
p1 = plot3(X, Y, Z);
hold on
p2 = plot3(X_pretie, Y_pretie, Z_pretie);

legend([p3 p1 p2],  'Actuated 6 deg','No control', 'Pretied 6 deg')
axis equal
title("NED Max Wind");
xlabel('N (m)'); ylabel('E (m)'); zlabel('D (m)');
hold off

load drift_control_nowind.mat
load drift_nocontrol_nowind.mat
load drift_pretie_nowind.mat

figure();
p3 = plot3(X_control, Y_control, Z_control);

hold on
p2 = plot3(X_pretie, Y_pretie, Z_pretie);
p1 = plot3(X, Y, Z);
legend([p3 p1 p2], 'Actuated 6 deg', 'No control', 'Pretied 6 deg')
axis equal
title("NED No Wind");
xlabel('N (m)'); ylabel('E (m)'); zlabel('D (m)');
hold off