function [x,y,err_max] = ERK(fid, a, b, h, F, eta, ARK, bRK, cRK, sol)
  %  ERK.m : Ejecuta el algoritmo ERK (Explicit Runge-Kutta) de s etapas
  %            para resolver una solución por aproximación discreta 
  %            del problema del valor inicial de una EDO en un intervalo [a,b].
  %
  %   ARGUMENTOS DE ENTRADA:
  %    -> fid : Control de fichero donde escribir resultados
  %    -> a : Valor inicial (x_0 = a)
  %    -> b : Extremo del intervalo a estudiar
  %    -> h : Paso entre iterantes
  %    -> F : Función que define la derivada (dydx(x) = F(x,y(x)))
  %    -> eta : Condición inicial (y(a) = eta)
  %    -> ARK : Matriz de coeficientes aij del tablero de Butcher del método.
  %    -> bRK : Vector de coeficientes bi del tablero de Butcher del método.
  %    -> cRK : Vector de coeficientes ci del tablero de Butcher del método.
  %    -> sol : Vector de funciones que son solución exacta del problema
  %
  %   SALIDA DE LA FUNCIÓN:
  %    -> x : Iterante final
  %    -> y : Solución numérica en el iterante final
  %    -> err_max : Error máximo cometido durante las iteraciones
  %
  
  N = floor((b - a) / h);
  m = length(eta);
  s = length(bRK);
  kn = zeros(m,s);
  
  x = a;
  y = eta;
  
  if nargin == 10
    err = y - sol(x);
    escribe_cabecera(fid, x, y, err)
    escribe_paso(fid, 0, x, y, err)
    err_max = err;
  else
    escribe_cabecera(fid, x, y)
    escribe_paso(fid, 0, x, y)
  endif
  
  for k=1:N     %bucle iter-total
    for i=1:s
      aux = y;
      for j=1:(i-1) %bucle kn's
        aux =  aux + h* ARK(i,j) * kn(:,j);
      endfor
      kn(:,i) = F(x + cRK(i)*h, aux);
    endfor
    for i=1:s   %bucle y_{n+1}
      y = y + h * bRK(i)*kn(:,i);
    endfor
    x = a + k * h;
    
    if nargin == 10
      err = y - sol(x);
      escribe_paso(fid, k, x, y, err)
      for i=1:length(err)
        err_max(i) = max(abs(err_max(i)), abs(err(i)));
      endfor
    else
      escribe_paso(fid, k, x, y)
    endif
  endfor
  
  fprintf(fid, '\n')
  fprintf(fid, 'Último instante iterado: x = \n %10.5e', x)
  fprintf(fid, '\n')
  fprintf(fid, 'Último iterante obtenido: y = ')
  fprintf(fid, '\n %10.5e', y)
  fprintf(fid, '\n')
  fprintf(fid, 'Error máximo cometido: err_max = ')
  fprintf(fid, '\n %10.5e', err_max)
  fprintf(fid, '\n')
  
  return
  