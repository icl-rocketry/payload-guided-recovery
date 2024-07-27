clear; close all
T = 0.05;
end_time = 1000;
num_of_steps = floor(end_time / T);

u=cell(1,num_of_steps);
x=cell(1,num_of_steps);

x{1} = [0; 0; 1000; 0; 0; 0; 7; 0; 2; 0; 0; 0];

for count = 1:num_of_steps
    u{count} = zeros(2,1);
    % if count * T > 14 && count * T < 14.5
    % u{count} = [0; 0.2];
    % end

    % if count * T > 3 && count * T < 4
    % u{count} = [0.75; 0];
    % end
end

disp = 0;

for i = 1:num_of_steps
    % x{step}(1:3)=[0;0;0];
    x{i+1}=RungeKutta4(@six_dof_parachute, x{i}, u{i}, T);
    if disp == 100
        disp = 0;
        T*i;
    else
        disp = disp + 1;
    end


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
    X(i) = x{i}(1);
    Y(i) = x{i}(2);
    Z(i) = x{i}(3);

    phi(i) = x{i}(4);
    theta(i) = x{i}(5);
    psi(i) = x{i}(6);

    U(i) = x{i}(7);
    V(i) = x{i}(8);
    W(i) = x{i}(9);

    P(i) = x{i}(10);
    Q(i) = x{i}(11);
    R(i) = x{i}(12);

    t(i) = i .* T;

    if x{1}(3) ~= 0
        if x{i+1}(3) <= 0
            break;
        end
    end
end




figure();
plot(t, X);

figure();
plot(t, U);

figure();
plot(t, Z);

% figure();
% plot3(X, Y, Z);

figure();
plot(t, W);