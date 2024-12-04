function [x,y,err_max] = EulerExp_draw(fid, a, b, h, F, eta, sol)
  % EulerExp_draw.m : Ejecuta el algoritmo de Euler explícito para calcular una
  %              aproximación discreta del problema del valor inicial de una EDO
  %              en un intervalo [a,b] dado. Además representa los resultados en
  %              una gráfica.
  %
  %   ARGUMENTOS EN EulerExp_draw
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
  fmt = ['bo'; 'ro'; 'go'];
  figure
  clf
  hold on
  
  if nargin == 7
    err = y - sol(x);
    escribe_cabecera(fid, x, y, err)
    escribe_paso(fid, 0, x, y, err)
    %err_max = norm(err, 2)
  else
    escribe_cabecera(fid, x, y)
    escribe_paso(fid, 0, x, y)
  endif
  for i=1:size(y)
      plot(x,y(i),fmt(i))
  endfor
  
  for k=1:N
    y = y + h * F(x,y);
    x = a + k * h;
    if nargin == 7
      err = y - sol(x);
      escribe_paso(fid, k, x, y, err)
      %err_max = max(err_max, err);
      for i=1:size(y)
        plot(x,y(i),fmt(i))
      endfor
    else
      escribe_paso(fid, k, x, y)
      for i=1:size(y)
        plot(x,y(i),fmt(i))
      endfor
    endif
  endfor
  
  hold off
  return
