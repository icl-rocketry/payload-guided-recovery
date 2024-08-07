function [x,u] = RungeKutta4(f,x,u,h)
% Fourth-order Runge-Kutta

k1 = f(x,u(:));
k2 = f(x+k1./2,u(:));
k3 = f(x+k2./2,u(:));
k4 = f(x+k3./2,u(:));
x(:,+1) = x(:) + h.*(k1/6 + (k2+k3)/3 + k4/6);
