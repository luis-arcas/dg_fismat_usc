function [x,y,err_max] = DIRK(fid, a, b, h, F, eta, ARK, bRK, cRK, sol)
  % DIRK.m : Ejecuta el algoritmo DIRK (Diagonally-Implicit Runge-Kutta) de s etapas
  %            para calcular una solución por aproximación discreta 
  %            del problema del valor inicial de una EDO definida en un intervalo [a,b] dado
  %
  %   ARGUMENTOS EN DIRK(fid, a, b, N, h, F, ARK, bRK, cRK, eta, sol)
  %    -> fid : Control de fichero donde escribir resultados
  %    -> a   : Valor inicial (x_0 = a)
  %    -> b   : Extremo del intervalo a estudiar
  %    -> h   : Paso entre iterantes
  %    -> F   : Función que define la derivada (dydx(x) = F(x,y(x)))
  %    -> eta : Condición inicial (y(a) = eta)
  %    -> ARK : Matriz de coeficientes aij del tablero de Butcher del método.
  %    -> bRK : Vector de coeficientes bi del tablero de Butcher del método.
  %    -> cRK : Vector de coeficientes ci del tablero de Butcher del método.
  %    -> sol : Vector de funciones que son solución exacta del problema
  %
  %   SALIDA DE LA FUNCIÓN
  %    -> x       : Iterante final
  %    -> y       : Solución numérica en el iterante final
  %    -> err_max : Error máximo cometido durante las iteraciones
  %
  
  N = floor((b - a) / h);
  m = length(eta);
  s = length(bRK);
  kn = zeros(m,s);
 
  x = a;
  y = eta;
  
    %--PRECISIÓN PARA IterBroyden--%
    eps1 = 1.e-8;
    eps2 = 1.e-10;
    itmax = 50;
  
  if nargin == 10
    err = y - sol(x);
    escribe_cabecera(fid, x, y, err)
    escribe_paso(fid, 0, x, y, err)
    err_max = err;
  else
    escribe_cabecera(fid, x, y)
    escribe_paso(fid, 0, x, y)
  endif
  
  %Fi = z - F(x,y)
  %--SOLUCIÓN POR fsolve--%
    %Kn1 = fsolve(F1, z, optimset('Display','off','TolFun',1.e-8,'TolX',1.e-8));
  %--SOLUCIÓN POR IterBroyden--%
     % -> A0 = D(F1)(x_{n+1},y_{n+1}) = Id - h*dFdy(x_{n+1},y_{n+1})
     % -> Si h<<1, es razonable aproximar A0~Id, como hacemos a continuación:
     
  A0 = eye(length(eta));
  for k=1:N
    Faux = @(z) z - F(x + cRK(1)*h, y + h * ARK(1, 1) * z);
    [kn(:,1), ~, ~] = IterBroyden(fid, F(x,y), Faux, A0, eps1, eps2, itmax);
    y = y + h * bRK(1) * kn(:,1);
    for i=2:s
      aux = y;
      for j=1:(i-1)
        aux = aux + h * ARK(i,j) * kn(:,j);
      endfor 
      Faux = @(z) z - F(x + cRK(i)*h, aux + h * ARK(i,i) * z);
      [kn(:,i), ~, ~] = IterBroyden(fid, kn(:,i-1), Faux, A0, eps1, eps2, itmax);
      y = y + h * bRK(i) * kn(:,i);
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