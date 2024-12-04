function [x,y] = EulerImp(fid, a, eta, h, N, F, sol)
  % EulerExp.m : Ejecuta el algoritmo de Euler expl�cito para calcular una
  %              aproximaci�n discreta del problema del valor inicial de una EDO
  %              en un intervalo [a,b] dado
  %
  %   ARGUMENTOS EN EulerExp(a, N, h, F, eta)
  %    -> fid : Control de fichero donde escribir resultados
  %    -> a : Valor inicial (x_0 = a)
  %    -> eta : Condici�n inicial (y(a) = eta)
  %    -> h : Paso entre iterantes
  %    -> N : N�mero total de pasos
  %    -> F : Campo que define la EDO (dydx(x) = F(x,y(x)))
  %    -> sol : Vector de funciones que son soluci�n exacta del problema
  %
  %   SALIDA DE LA FUNCI�N
  %    -> x : Iterante final
  %    -> y : Soluci�n num�rica en el iterante final
  %    -> err_max : Error m�ximo cometido durante las iteraciones
  %
  
  x = a;
  y = eta;
  %--PRECISI�N PARA IterBroyden--%
    eps1 = 1.e-8;
    eps2 = 1.e-10;
    itmax = 50;
  
  
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
    x = a + k * h;
    FEI = @(z) z - y + h * F(x,z);
    %--SOLUCI�N POR fsolve--%
    %y = fsolve(FEI, y, optimset('Display','off','TolFun',1.e-8,'TolX',1.e-8));
    %--SOLUCI�N POR IterBroyden--%
     % -> A0 = D(FEI)(x_{n+1},y_{n+1}) = Id - h*dFdy(x_{n+1},y_{n+1})
     % -> Si h<<1, es razonable aproximar A0~Id, como hacemos a continuaci�n:
    A0 = eye(length(eta));
    [y, ~, ~] = IterBroyden(fid, y, FEI, A0, eps1, eps2, itmax);
    if nargin == 7
      err = y - sol(x);
      escribe_paso(fid, k, x, y, err)
      %err_max = max(err_max, err);
    else
      escribe_paso(fid, k, x, y)
    endif
  endfor
  fprintf(fid,'Tiempo final =  %10.3e \n', x)
  fprintf(fid,'Soluci�n num�rica =  %10.5e \n', y)
  %Implementaci�n alternativa
  %%FEI = (Z-y)1/h - F(x,z); A0 = 1/h - dFdy(x,y) (gradiente en segunda variable)
  return
