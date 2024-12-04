format long
clear all

fid = 1;
F = @(x) sin(x) - x^2;
A0 = -1.;
x0 = 0.5;
eps = 1.e-6;
eta = 1.e-6;
itmax = 30;

[x,A,k] = IterBroyden(fid, x0, F, A0, eps, eta, itmax)