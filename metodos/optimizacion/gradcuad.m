function [x,norm_grad,index] = gradcuad(fid, A, b, x0, eps, itmax, sol)

y = @(x) 0.5 * transpose(x) * A * x - transpose(b) * x;
x = x0;
g = A*x - b;

  %if nargin == 7
  %  err = y(x) - sol;
  %  escribe_cabecera(fid, x, y(x), err)
  %  escribe_paso(fid, 0, x, y(x), err)
  %  err_max = err;
  %else
  %  escribe_cabecera(fid, x, y(x))
  %  escribe_paso(fid, 0, x, y(x))
  %endif
  
  norm_grad = norm(g,2);
  if norm_grad <= eps
    fprintf(fid, 'Gradiente muy próximo a cero')
    index = 0;
    return
  endif

for k=1:itmax
  Ag = A * g;
  gAg = transpose(g) * Ag;
  aopt = norm_grad^2 / gAg;
  x = x - aopt*g;

  %if nargin == 7
  %    err = y(x) - sol;
  %    escribe_paso(fid, k, x, y(x), err)
  %    for i=1:length(err)
  %      err_max(i) = max(abs(err_max(i)), abs(err(i)));
  %    endfor
  %  else
  %    escribe_paso(fid, k, x, y(x))
  %endif
  
  %g = grad(x);
  g = g - aopt * Ag;
  norm_grad = norm(g,2);
  if norm_grad <= eps
    fprintf(fid, 'Gradiente muy próximo a cero')
    index = k;
    return
  endif
endfor
    fprintf(fid, 'Iteraciones máximas alcanzadas sin convergencia')
    index = -1;
return