clear all
format long
%fid = 1;
fid = fopen('SalidaAB3.txt', 'w');

###################################
###MLM ADAMS-BASHFORT DE 3 PASOS###
###################################


#Ejercicio 1 - Caso escalar

fprintf(fid, '########################## \n')
fprintf(fid, 'EJERCICIO 1 - CASO ESCALAR \n')
fprintf(fid, '########################## \n')
fprintf(fid, '\n')

a = 2.;
b = 3.;

F = @(x,y) 1. + (x-y)^2;
exact = @(x) x + (1 - x)^-1;
eta = 1.;
       
fprintf(fid, 'Longitud de paso h = %10.3e', 1.e-1)
fprintf(fid, '\n')
[x,y,err_max] = AB3(fid, a, b, 0.1, F, eta, exact)
fprintf(fid, '\n')


fprintf(fid, 'Longitud de paso h = %10.3e', 1.e-2)
fprintf(fid, '\n')
[x,y,err_max] = AB3(fid, a, b, 0.01, F, eta, exact)
fprintf(fid, '\n')


fprintf(fid, 'Longitud de paso h = %10.3e', 1.e-3)
fprintf(fid, '\n')
[x,y,err_max] = AB3(fid, a, b, 0.001, F, eta, exact)
fprintf(fid, '\n')


#Ejercicio 2 - Caso vectorial

fprintf(fid, '############################ \n')
fprintf(fid, 'EJERCICIO 2 - CASO VECTORIAL \n')
fprintf(fid, '############################ \n')
a = 0.;
b = 2.4;


F = @(x,y) [y(2); -y(1) + cos(3.*x)];
exact = @(x) [(9.*cos(x) - cos(3.*x)) / 8.; 
              3. * (-3.* sin(x) + sin(3.*x)) / 8.];
eta = [1.; 0.];

h = 8.e-2;
fprintf(fid, 'Longitud de paso h = %10.3e', h)
fprintf(fid, '\n')
[x,y,err_max] = AB3(fid, a, b, h, F, eta, exact)
fprintf(fid, '\n')



#Ejercicio 3 - Sistema de Lorenz

%fprintf(fid, '############################')
%fprintf(fid, 'EJERCICIO 3 - SISTEMA LORENZ')
%fprintf(fid, '############################')

#HACER!!!!!!

%lambda = -1.e4;
%x_star = log(2) / abs(lambda)
%a = 0.;
%b = 10. * x_star;
%h = x_star / 10.;

%F = @(x,y) cos(x) + lambda * (y - sin(x));
%exact = @(x) exp(lambda * x) + sin(x);
%eta = 1.;

%fprintf(fid, '\n')
%[x,y] = EulerExp(fid, a, b, h, F, eta, exact)

%fprintf(fid, 'Introduciendo nuevo h, h_2: \n')
%h2 = 1.e-4;
%[x2,y2] = EulerExp(fid, x, 2*pi, h2, F, y, exact)

%Con h2 requerido para la estabilidad numérica, se necesitan más de 60,000
%iteraciones! Para nada viable. Necesitamos otro método. Nos pasamos a EulerImp!

%########## IMPORTANTE SI FID != 1 ############
fclose(fid);



