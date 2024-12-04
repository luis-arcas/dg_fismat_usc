clear
clc
format long
#fid = 1;
fid = fopen('SalidaEulerImp.txt', 'w');

#####################
###EULER IMPLÍCITO###
#####################

#Ejercicio 1 - Caso escalar

%fprintf(fid, '##########################')
%fprintf(fid, 'EJERCICIO 1 - CASO ESCALAR')
%fprintf(fid, '##########################')
a = 2.;
b = 3.;
h = 0.01;
N = floor((b - a) / h);

F = @(x,y) 1. + (x-y)^2;
exact = @(x) x + 1. / (1. - x);
eta = 1.;

fprintf(fid, '\n')
[x,y] = EulerImp(fid, a, eta, h, N, F, exact);
fclose(fid);
