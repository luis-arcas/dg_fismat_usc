format long
clear all

F = @(x) sin(x) - x^2;
dFdx = @(x) cos(x) - 2*x;
x0 = 1.;
eps = 1.e-6;
eta = 1.e-6;
itmax = 30;

[x, sol, dif] = IterNewton(x0, F, dFdx, eps, eta, itmax)