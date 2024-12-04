clear
clc
format long

fid = fopen('SalidaMainGrad.txt', 'w');
n = 5;
A = diag(4*ones(1,n),0) + ...
    diag(ones(1,n-1),1) + ...
    diag(ones(1,n-1),-1);

b = zeros(n,1);
b(1) = 3;
for i=2:n-1
  k = mod(i,2);
  b(i) = (-1)^k * 2;
endfor
b(n) = (-1)^(n+1) * 3;
x0 = zeros(n,1);

eps = 1.e-6;
itmax = 1000;
%sol = ;

[x,norm_grad,index] = gradcuad(fid, A, b, x0, eps, itmax)
fclose(fid);