open("SimulationGuidedRecoveryV2.slx")
%open("SimulationGuidedRecoveryV2_with_drag.slx")
out = sim("SimulationGuidedRecoveryV2.slx");

X = out.X.data;
Y = out.Y.data;
Z = out.Z.data;

figure(1)
clf;
index = find(Z == 0, 1);
X = X(1:index);
Y = Y(1:index);
Z = Z(1:index);

plot3(X,Y,Z,'LineWidth',3)
hold on
scatter3(0,0,0)


hold on
Vvec = out.Vvec.data();
Vx = Vvec((1:100:index),1);
Vy = Vvec((1:100:index),2);
Vz = Vvec((1:100:index),3);
quiver3(X(1:100:end),Y(1:100:end),Z(1:100:end),Vx,Vy,Vz, 'AutoScale', 'on')
hold on
V1 = cos(out.ErrorHeading.data);
V2 = sin(out.ErrorHeading.data);
V3 = zeros(size(V2,1),1);
%quiver3(X,Y,Z,V1,V2,V3, 'AutoScale', 'on')
hold on
V11 = out.VecToTar.data(:,1);
V22 = out.VecToTar.data(:,2);
%quiver3(X,Y,Z,V11,V22,V3, 'AutoScale', 'on')

axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');

xlabel("x")
ylabel("y")
zlabel("z")

xlim([min(X) max(X)])
ylim([min(Y) max(Y)])

figure(2)
clf;
VMag = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
VMag = VMag(1:index);
plot((1:1:index),VMag)