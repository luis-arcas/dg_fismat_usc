function [x,y,err_max] = EulerExp(fid, a, b, h, F, eta, sol)
  % EulerExp.m : Ejecuta el algoritmo de Euler expl�cito para calcular una
  %              aproximaci�n discreta del problema del valor inicial de una EDO
  %              en un intervalo [a,b] dado
  %
  %   ARGUMENTOS EN EulerExp(a, N, h, F, eta)
  %    -> fid : Control de fichero donde escribir resultados
  %    -> a : Valor inicial (x_0 = a)
  %    -> b : Extremo del intervalo a estudiar
  %    -> h : Paso entre iterantes
  %    -> F : Funci�n que define la derivada (dydx(x) = F(x,y(x)))
  %    -> eta : Condici�n inicial (y(a) = eta)
  %    -> sol : Vector de funciones que son soluci�n exacta del problema
  %
  %   SALIDA DE LA FUNCI�N
  %    -> x : Iterante final
  %    -> y : Soluci�n num�rica en el iterante final
  %    -> err_max : Error m�ximo cometido durante las iteraciones
  %
  
  N = floor((b - a) / h);
  x = a;
  y = eta;
  
  
  if nargin == 7
    err = y - sol(x);
    escribe_cabecera(fid, x, y, err)
    escribe_paso(fid, 0, x, y, err)
    %err_max = norm(err, 2)
  else
    escribe_cabecera(fid, x, y)
    escribe_paso(fid, 0, x, y)
  endif
  
  for k=1:N
    y = y + h * F(x,y);
    x = a + k * h;
    if nargin == 7
      err = y - sol(x);
      escribe_paso(fid, k, x, y, err)
      %err_max = max(err_max, err);
    else
      escribe_paso(fid, k, x, y)
    endif
    
  endfor
  
  return
