function [x,u] = RungeKutta4(f,x,u,h)
% Fourth-order Runge-Kutta

k1 = f(x,u(:));
k2 = f(x+k1./2*h,u(:));
k3 = f(x+k2./2*h,u(:));
k4 = f(x+k3./2*h,u(:));
x(:,+1) = x(:) + h.*(k1/6 + (k2+k3)/3 + k4/6);

% k_1 = func(t_0    , y_0          );
% k_2 = func(t_0+h/2, y_0+(h/2)*k_1);
% k_3 = func(t_0+h/2, y_0+(h/2)*k_2);
% k_4 = func(t_0+h  , y_0+    h*k_3);
% y = y_0 + (h/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);