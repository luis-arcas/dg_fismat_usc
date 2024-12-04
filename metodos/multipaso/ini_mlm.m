function [y,f,err_max] = ini_mlm(fid, a, N, h, F, eta, sol)
  % ini_mlm.m : Ejecuta N iteraciones del algoritmo de Heun (RungeKutta de orden 2)
  %             para inicializar un método lineal multipaso (MLM) de orden de convergencia 2.
  %
  %   ARGUMENTOS EN EulerExp(a, N, h, F, eta)
  %    -> fid : Control de fichero donde escribir resultados
  %    -> a : Valor inicial (x_0 = a)
  %    -> N : Número de iteraciones requeridas
  %    -> h : Paso entre iterantes
  %    -> F : Función que define la derivada (dydx(x) = F(x,y(x)))
  %    -> eta : Condición inicial (y(a) = eta)
  %    -> sol : Vector de funciones que son solución exacta del problema
  %
  %   SALIDA DE LA FUNCIÓN
  %    -> y : Solución numérica en el iterantes calculados
  %    -> f : Vector de evaluaciones de F en los iterantes calculados
  %    -> err_max : Error máximo cometido durante las iteraciones
  %
  
  x = a;
  y(:,1) = eta;
  f(:,1) = F(x,y(:,1));
  
  if nargin == 7
    err = y(:,1) - sol(x);
    escribe_cabecera(fid, x, y(:,1), err)
    escribe_paso(fid, 0, x, y(:,1), err)
    err_max = err;
  else
    escribe_cabecera(fid, x, y(:,1))
    escribe_paso(fid, 0, x, y(:,1))
  endif
  
  for k=2:N
    Kn1 = f(:,k-1);
    Kn2 = F(x + h, y(:,k-1) + h * Kn1);
    y(:,k) = y(:,k-1) + h / 2 * (Kn1 + Kn2);
    x = a + (k-1) * h;
    f(:,k) = F(x,y(:,k));
    if nargin == 7
      err = y(:,k) - sol(x);
      escribe_paso(fid, k, x, y(:,k), err)
      for i=1:length(err)
        err_max(i) = max(abs(err_max(i)), abs(err(i)));
      endfor
    else
      escribe_paso(fid, k, x, y(:,k))
    endif
   
  endfor
  #PUNTO ESTRELLA (.*) = PRODUCTO COMPONENTE A COMPONENTE DE MATRICES
  return
