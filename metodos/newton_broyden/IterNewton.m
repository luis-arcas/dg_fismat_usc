function [x, sol, dif] = IterNewton(x0, F, dFdx, eps, eta, itmax)
  x = x0;
  sol = fsolve(F,x0);

  for k=1:itmax
    A = dFdx(x);
    b = F(x);
    delta = - A \ b;
    x = x + delta;
    if (norm(delta,2)<eps)
      dif = norm(x-sol,2);
      fprintf(1,"Parada test (iterantes muy próximos).");
        fprintf(1, '\n');
      return
    elseif (norm(F(x),2)<eta)
      dif = norm(x-sol,2);
      fprintf(1,"Parada test (solución muy próxima).");
        fprintf(1, '\n');
      return
    endif
  endfor
  dif = norm(x-sol,2);
  fprintf(1,"Parada test (n.º de iter. máximas alcanzado).");
    fprintf(1, '\n');
 return

