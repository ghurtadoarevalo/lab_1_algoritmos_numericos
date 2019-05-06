function [sol,errors,iterations] = gauss_seidel(A, b)
    % sol es vector solución
    % err es el error final
    % it es el numero de iteraciones final
    % A matriz del sistema
    % b vector del sistema
    % maxiter es numero máximo de iteraciones
    % epsilon es la cota del error
    
    n=length(b);
    sol=zeros(1,n)';

    error = norm(A*sol - b);

    tolerance = 0.00005;
    iterations = 0;
    maxiter = 100;

    errors = [];
    errors = [errors, error];

    while error > tolerance
        sol_old = sol;  
        for i=1:n
            S=0;
            for j=1:i-1
                S=S+A(i,j)*sol(j);
            end

            for j=i+1:n
                S=S+A(i,j)*sol_old(j);
            end

            sol(i)=(1/A(i,i))*(b(i)-S);
        end

        iterations=iterations+1;
        error=norm(sol_old-sol);
        errors = [errors, error];
    end
end