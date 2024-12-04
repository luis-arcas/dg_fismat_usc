function [x,y,err_max] = ERK4(fid, a, b, h, F, eta, sol)
  % ERK4.m : Ejecuta el algoritmo de RungeKutta de orden 4 en k's para calcular una
  %              aproximación discreta del problema del valor inicial de una EDO
  %              en un intervalo [a,b] dado
  %
  %   ARGUMENTOS EN EulerExp(a, N, h, F, eta)
  %    -> fid : Control de fichero donde escribir resultados
  %    -> a : Valor inicial (x_0 = a)
  %    -> b : Extremo del intervalo a estudiar
  %    -> h : Paso entre iterantes
  %    -> F : Función que define la derivada (dydx(x) = F(x,y(x)))
  %    -> eta : Condición inicial (y(a) = eta)
  %    -> sol : Vector de funciones que son solución exacta del problema
  %
  %   SALIDA DE LA FUNCIÓN
  %    -> x : Iterante final
  %    -> y : Solución numérica en el iterante final
  %    -> err_max : Error máximo cometido durante las iteraciones
  %
  
  N = floor((b - a) / h);
  x = a;
  y = eta;
  
  
  if nargin == 7
    err = y - sol(x);
    escribe_cabecera(fid, x, y, err)
    escribe_paso(fid, 0, x, y, err)
    err_max = err;
  else
    escribe_cabecera(fid, x, y)
    escribe_paso(fid, 0, x, y)
  endif
  
  for k=1:N
    Kn1 = F(x,y);
    Kn2 = F(x + h / 2, y + h * Kn1 / 2);
    Kn3 = F(x + h / 2, y + h * Kn2 / 2);
    Kn4 = F(x + h, y + h * Kn3);
    y = y + h / 6 * (Kn1 + 2*Kn2 + 2*Kn3 + Kn4);
    x = a + k * h;
    if nargin == 7
      err = y - sol(x);
      escribe_paso(fid, k, x, y, err)
      for i=1:length(err)
        err_max(i) = max(abs(err_max(i)), abs(err(i)));
      endfor
    else
      escribe_paso(fid, k, x, y)
    endif
  endfor
  #PUNTO ESTRELLA (.*) = PRODUCTO COMPONENTE A COMPONENTE DE MATRICES
  return
