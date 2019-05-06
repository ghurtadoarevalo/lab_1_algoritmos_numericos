function [sol,errors,iterations] = gauss_jacobi(A, b)  
    %Si A es una matriz de diagonal estrictamente dominante por filas,
    %El método converge.

    maxiter = 100;
    epsilon = 0.000000005;
    
    % sol es vector solución
    % error es el error final
    % it es el numero de iteraciones final
    % A matriz del sistema
    % b vector del sistema
    % maxiter es numero máximo de iteraciones
    % epsilon es la cota del error
    
    n=length(b);
    iterations=0;
    sol_x=zeros(1,n)';
    sol=zeros(1,n);

    error = norm(A*sol_x - b);
    errors = [];
    errors = [errors, error];

    while error>epsilon
        for i=1:n
           S=0;
           for j=1:n
               if i~=j
                S=S+A(i,j)*sol_x(j);
               end
           end
           sol(i)=(b(i)-S)/A(i,i);
        end
        iterations=iterations+1;
        error=norm(sol_x-sol);
        errors = [errors, error];
        sol_x=sol;
    end
end