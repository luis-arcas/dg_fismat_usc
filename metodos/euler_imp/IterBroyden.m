function [x,A,k] = IterBroyden(fid, x0, F, A0, eps, eta, itmax)
  %IterBroyden: Función que permite calcular por iteración la raíz de una función
  % a partir de dos iterantes iniciales y una matriz de paso inicial que aproxima
  % a la jacobiana de la función.
  %
  %
  %
  %
  
  x = x0;
  A = A0;
  k = 0;
  
   for k=1:itmax
     
     b = -F(x);
     s = A\b;
     x = x + s;
     y = F(x) - F(x - s);
     mat = y - A * s;
     mat = mat * transpose(s);
     mat = mat / norm(s,2)^2;
     A = A + mat;
     if (norm(F(x),2)<eta)
       fprintf(fid, "Parada test (solución muy próxima) \n")
       return
     endif
     if (norm(s,2)<eps)
       fprintf(fid, "Parada test (iterantes muy próximos) \n")
       return
     endif
   endfor
  fprintf(fid, "Parada test (itmax alcanzado) \n")
  return