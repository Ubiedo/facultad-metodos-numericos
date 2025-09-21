% Federico Javier Gonzalez Ubiedo -------------------------------------------- %
% Parte c, Ejercicio 12, Practico 2 ------------------------------------------ %
A = [1 1 1 1; -.2 .6 -.3 -.1; -.2 -.2 .6 -.2; -.1 -.1 -.1 .4]; % (I - P^t) --- %
b = [1; 0; 0; 0];
x = A \ b; % da un vector x con la solucion computacional al sistema Ax = b -- %
for i = 1:4 % bucle para ver la solucion ------------------------------------- %
  fprintf('x%d = %f\n', i, x(i));
endfor;
