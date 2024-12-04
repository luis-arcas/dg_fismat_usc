% PREDICTOR-CORRECTOR!!

% P      : Predigo un iterante inicial y_n+1 = y_n+1(AB3) para AM2.
% (EC)^L : Evalúo y compruebo si es un buen iterante. Podemos repetir este proceso
%          tantas veces como sea necesario (L comprobaciones).
% E      : Una vez obtengo un buen iterante, evalúo y ejecuto AM2.
% ESTRUCTURA: P (EC)^L E

yold = y
x = a + h*k
y = yold + AM3!!

for 