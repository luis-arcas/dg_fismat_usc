function [x,y,err_max] = AM2(fid, a, b, h, F, eta, sol)
  % AB3.m : Ejecuta el algoritmo de Adam Bashfort de k=3 pasos para resolver el
  %         el problema del valor inicial y' = F(x,y) en un intervalo [a,b] dado. 
  %
  %   ARGUMENTOS EN EulerExp(a, N, h, F, eta)
  %    -> fid : Control de fichero donde escribir resultados
  %    -> a : Valor inicial (x_0 = a)
  %    -> b : Último iterante a calcular.
  %    -> h : Paso entre iterantes
  %    -> F : Función que define la derivada (dydx(x) = F(x,y(x)))
  %    -> eta : Condición inicial (y(a) = eta)
  %    -> sol : Vector de funciones que son solución exacta del problema
  %
  %   SALIDA DE LA FUNCIÓN
  %    -> x : Último iterante calculado.
  %    -> y : Solución numérica en el último iterante.
  %    -> err_max : Error máximo cometido durante las iteraciones
  %
  
  N = floor((b - a) / h);
  kAM = 2;
  [y_ini, f, err_max] = ini_mlm(fid, a, kAM, h, F, eta, sol);
  y = y_ini(kAM);
  %options = ...;
  
  for k=kAM+1:N
    FIMP = @(z) z - y - h/12 * (5*F(x+h*k, z) + 8*f(:,2) - f(:,1));
    [~,y,~] = IterBroyden(...)
    x = a + h * k;
    f(:,1) = f(:,2);
    f(:,2) = f(:,3);
    f(:,3) = F(x,y);
    
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