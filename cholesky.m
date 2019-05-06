function [sol,error,costo]=cholesky(A, b)

    % sol es vector soluci�n
    % error es el error final
    % A matriz del sistema
    % b vector del sistema

    [n,m]=size(A);
    costo = 0;
    for k = 1 : n
        A(k,k) = sqrt(A(k,k));
        A(k+1:n,k) = A(k+1:n,k)/A(k,k);
        costo= costo+2;
        for j = k + 1 : n
            A(j:n,j) = A(j:n,j) - A(j,k)*A(j:n,k);
            costo = costo + 2;
        end
    end
    
    L = tril(A);
    z = inv(L)*b;
    sol = inv(L')*z;
    error = norm(eye(n)-inv(L*L')*A);
    costo = costo + 1;
end