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
eta = [0.; 0.];  % y(a).
f   = @(x,y) [y(2); 2 * x * sin(x) - y(1)]; % Funcion que define la EDO.
sol = @(x) [x / 2 * sin(x) - x^2 / 2 * cos(x);
            sin(x) * (1 + x^2) / 2 - 1/2 * x * cos(x)];   % Solucion del PVI.

% Parametros para el integrador.
h   = 0.01;        % Parametro de discretizacion.
N   = floor((b - a) / h);        % Numero total de pasos.
fid = fopen('SalidaEj1.txt','w'); % fid para salidas.

% Llamada al Runge-Kutta. 
 [x,y] = ERK3(fid, f, a, eta, h, N, sol);

% Cierre del fichero de salida.
fclose(fid);





