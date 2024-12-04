clear all
format long
fid = 1;
%fid = fopen('SalidaEuler1.txt', 'w');

#####################
###EULER EXPL�CITO###
#####################

#Ejercicio 1 - Caso escalar

%fprintf(fid, '##########################')
%fprintf(fid, 'EJERCICIO 1 - CASO ESCALAR')
%fprintf(fid, '##########################')
%a = 2.;
%b = 3.;
%h = 1.e-3;

%F = @(x,y) 1. + (x-y)^2;
%exact = @(x) x + (1 - x)^-1;
%eta = 1.;

%fprintf(fid, '\n')
%[x,y,err_max] = EulerExp(fid, a, b, h, F, eta, exact)


#Ejercicio 2 - Caso vectorial

%fprintf(fid, '###########################')
%fprintf(fid, 'EJERCICIO 2 - CASO VECTORIAL')
%fprintf(fid, '############################')
%a = 0.;
%b = 1.;
%h = 8.e-2;

%F = @(x,y) [y(2); -y(1) + cos(3.*x)];
%exact = @(x) [(9.*cos(x) - cos(3.*x)) / 8.; 
%              3. * (-3.* sin(x) + sin(3.*x)) / 8.];
%eta = [1.; 0.];

%fprintf(fid, '\n')
%[x,y] = EulerExp(fid, a, b, h, F, eta, exact)


#Ejercicio 3 - Problema de estabilidad num�rica

%fprintf(fid, '##################################')
%fprintf(fid, 'EJERCICIO 3 - ESTABILIDAD NUM�RICA')
%fprintf(fid, '##################################')
lambda = -1.e4;
x_star = log(2) / abs(lambda)
a = 0.;
b = 10. * x_star;
h = x_star / 10.;

F = @(x,y) cos(x) + lambda * (y - sin(x));
exact = @(x) exp(lambda * x) + sin(x);
eta = 1.;

fprintf(fid, '\n')
[x,y] = EulerExp(fid, a, b, h, F, eta, exact)

fprintf(fid, 'Introduciendo nuevo h, h_2: \n')
h2 = 1.e-4;
[x2,y2] = EulerExp(fid, x, 2*pi, h2, F, y, exact)

%Con h2 requerido para la estabilidad num�rica, se necesitan m�s de 60,000
%iteraciones! Para nada viable. Necesitamos otro m�todo. Nos pasamos a EulerImp!

%########## IMPORTANTE SI FID != 1 ############
%fclose(fid);



