% se define la funcion f del ejercicio 10 ------------------------------------------------- %
f = @(x) log(x);
puntos = [0.1; 1; 2; 2.9];
% parte a --------------------------------------------------------------------------------- %
% interpolar f con un polinomio c´ubico por los nodos de abscisas x = 0.1, 1, 2, 2.9        %
% evaluar la interpolante en x = 1.5 y determinar el error de interpolaci´on en dicho punto %
A = [1, 0.1, 0.01, 0.001;
     1, 1, 1, 1;
     1, 2, 4, 8;
     1, 2.9, 8.41, 24.389];
b = [f(0.1); 0; f(2); f(2.9)];

% Ax = b                                                                                    %
function x = escal_gauss_pivoteo(A, b)
    [m, n] = size(A);
    Ab = [A, b];
    % escalonamiento gausiano con pivoteo parcial
    for k = 1:n-1
        [~, pivot] = max(abs(Ab(k:m, k)));
        pivot = pivot + k - 1;
        % intercambiar fila
        if pivot ~= k
            Ab([k, pivot], :) = Ab([pivot, k], :);
        endif
        for i = k+1:m
            factor = Ab(i, k) / Ab(k, k);
            Ab(i, k:n+1) = Ab(i, k:n+1) - factor * Ab(k, k:n+1);
        endfor
    endfor
    % sustitcion hacia atras
    x = zeros(n, 1);
    x(n) = Ab(n, n+1) / Ab(n, n);
    for i = n-1:-1:1
        x(i) = (Ab(i, n+1) - Ab(i, i+1:n) * x(i+1:n)) / Ab(i, i);
    endfor
endfunction

egp = escal_gauss_pivoteo(A, b);
% funcion interpolada de parte a                                                            %
pol_3 = @(x) egp(4) * x^3 + egp(3) * x^2 + egp(2) * x + egp(1);
% su valor en 1.5                                                                           %
disp('La interpolante en 1.5 es: ');
disp(pol_3(1.5));
% error de interpolacion en dicho punto                                                    %
error_int = f(1.5) - pol_3(1.5);
disp('El error de interpolacion en 1.5 es: ');
disp(error_int);
% parte b --------------------------------------------------------------------------------- %
% interpolar f mediante la interpolante de Hermite por x = 1, 2. Evaluar la interpolante en x = 1,5
% y determinar el error de interpolaci´on en dicho punto.
y = b;
d = [10; 1; 0.5; 10/29];
% h = 2 - 1 = 1, por lo que no lo tengo en cuenta en multiplicaciones y diviciones
pol_h = @(x) (((3*(x-1)^2)-2*(x-1)^3)*f(2)) + ((((x-1)^2)*(x-2))/2) + ((x-1)*((x-2)^2));                                                 %
% su valor en 1.5                                                                           %
disp('La interpolante en 1.5 es: ');
disp(pol_3(1.5));
% su error de interpolacion                                                                 %
error_her = ((((0.5)^2)*((-0.5)^2))/24)*((-6)/(1.5^4))
disp('El error de interpolacion de hermite en 1.5 es: ');
disp(error_her);
% parte c --------------------------------------------------------------------------------- %
x_range = linspace(1, 2, 200);
f_vals = arrayfun(f, x_range);
pol_3_vals = arrayfun(pol_3, x_range);
pol_h_vals = arrayfun(pol_h, x_range);
% graficar funciones                                                                  %
figure;
plot(x_range, f_vals, 'g', 'LineWidth', 1);
hold on;
plot(x_range, pol_3_vals, 'b', 'LineWidth', 1);
hold on;
plot(x_range, pol_h_vals, 'r', 'LineWidth', 1);

xlabel('x');
ylabel('f(x), pol_3(x), pol_h(x)');
title('Comparación entre f(x) y su Aproximación Polinómica');
legend('f(x) : log(x)', 'a : pol\_3(x)', 'b : pol\_h(x)');
grid on;
hold off;
