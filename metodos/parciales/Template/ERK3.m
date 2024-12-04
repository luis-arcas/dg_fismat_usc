function [x,y] = ERK3(fid,f,a,eta,h,N,ysol)

%
% ERK3.m
%
% Funcion destinada a integrar numericamente un pvi del tipo
%
%    dy(x)/dx = f(x,y(x)),    x in [a,b]
%        y(a) = eta,
%
% mediante un integrador temporal Runge-Kutta explicito de tres etapas.
% 
% Argumentos de entrada:
%    fid  : Identificador de fichero para escribir los resultados.
%    f    : Funcion que define la EDO.
%    a    : Tiempo inicial.
%    eta  : Valor de la solucion en el instante inicial.
%    h    : Paso de discretizacion.
%    N    : Numero total de pasos para llegar al instante final.
%    ysol : Solucion exacta del PVI.
%
% Argumentos de salida:
%    x    : Tiempo final de integracion.
%    y    : Solucion aproximada en el instante final.
%

% Inicializacion del integrador.
y   = eta;
x   = a;

% Calculo del error cometido y salida usando fid.
err = ysol(x)-y;
escribe_cabecera(fid,x,y,err);
escribe_paso(fid,0,x,y,err);

% Bucle en pasos.
for n=1:N

% Actualizacion de la solucion aproximada y del tiempo.
    Kn1 = f(x,y);
    Kn2 = f(x + h, y + h * Kn1);
    Kn3 = f(x + 2*h/3, y + 2*h/9 * (2*Kn1 + Kn2));
    
    y = y + h/4 * (Kn1 + 3*Kn3);
    x = a + n * h;
    
% Calculo del error cometido y salida usando fid.
    err = ysol(x) - y;
    escribe_paso(fid,n,x,y,err)
    
end

end