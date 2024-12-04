function [x,y,err_max] = DIRK2(fid, a, b, h, F, eta, sol)
  % DIRK2.m : Ejecuta el algoritmo DIRK (Diagonally-Implicit Runge-Kutta) de orden 2 en k's
  %            para calcular una solución por aproximación discreta 
  %            del problema del valor inicial de una EDO definida en un intervalo [a,b] dado
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
    %--PRECISIÓN PARA IterBroyden--%
    eps1 = 1.e-8;
    eps2 = 1.e-10;
    itmax = 50;
  
  
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
    %--SOLUCIÓN POR fsolve--%
    %Kn1 = fsolve(F1, z, optimset('Display','off','TolFun',1.e-8,'TolX',1.e-8));
    %--SOLUCIÓN POR IterBroyden--%
     % -> A0 = D(F1)(x_{n+1},y_{n+1}) = Id - h*dFdy(x_{n+1},y_{n+1})
     % -> Si h<<1, es razonable aproximar A0~Id, como hacemos a continuación:
    F1 = @(z) z - F(x + h/4, y + h/4 * z);
    %Kn1 = fsolve(F1, F(x,y), optimset('Display','off','TolFun',1.e-8,'TolX',1.e-8));
    A0 = eye(length(eta));
    [Kn1, ~, ~] = IterBroyden(fid, y, F1, A0, eps1, eps2, itmax);
    F2 = @(z) z - F(x + 3*h/4, y + h/2 * Kn1 + h/4 * z);
    %Kn2 = fsolve(F2, F(x,y+Kn1), optimset('Display','off','TolFun',1.e-8,'TolX',1.e-8));
    [Kn2, ~, ~] = IterBroyden(fid, y, F2, A0, eps1, eps2, itmax);
    y = y + h / 2 * (Kn1 + Kn2);
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