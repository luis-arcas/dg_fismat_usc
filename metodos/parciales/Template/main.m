% Programa principal destinado a resolver numericamente 
% un PVI del tipo 
%
%    dy(x)/dx = f(x,y(x)),    x in [a,b]
%        y(a) = eta,
%
% mediante un integrador temporal de tipo Runge-Kutta.

% Definicion del PVI.
a   = 0.;        % Tiempo inicial.
b   = 2.;        % Tiempo final.
eta = [6.; 1.; -1.];        % y(a).
f   = @(x,y) [y(2); 
              3*x^2 - x - y(2); 
              6*x - 1 - y(3)]; % Funcion que define la EDO.
sol = @(x) [6*exp(-x) + 7*x - 3.5 * x^2 + x^3; 
            -6*exp(-x) + 7. - 7*x + 3*x^2;
            6*exp(-x) - 7. + 6*x];   % Solucion del PVI.

% Parametros para el integrador.
h   = 1.e-3;        % Parametro de discretizacion.
N   = floor((b-a)/h);        % Numero total de pasos.
fid = fopen('SalidaEj1.txt','w'); % fid para salidas.

% Llamada al Runge-Kutta. 
[x,y] = ERK3(fid, f, a, eta, h, N, sol);

% Cierre del fichero de salida.
fclose(fid);





