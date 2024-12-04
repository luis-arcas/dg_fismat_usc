#Ejercicio 3 - MODELO SIR

fid = 1;

fprintf(fid, '########################')
fprintf(fid, 'EJERCICIO 3 - MODELO SIR')
fprintf(fid, '########################')
a = 0.;
b = 10.;
h = 1.e-1;

beta = 0.5;
gamma = 0.8;

F = @(y) [-beta * y(1) * y(2); 
          beta * y(1) * y(2) - gamma * y(2);
          gamma * y(2)];

total = 976.;
inf_ini = total^(-1);
rec_ini = 0.;
eta = [1. - inf_ini - rec_ini; inf_ini; rec_ini]

fprintf(fid, '\n')
[x,y] = EulerExp_draw(fid, a, b, h, F, eta)
fclose(fid);