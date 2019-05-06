function [sol,error,contador] = householder(A,b)

    % sol es vector solución
    % error es el error final
    % A matriz del sistema
    % b vector del sistema

    contador = 0;
    [m,n] = size(A);
    Q = eye(m);      
    R = A;
    for j = 1:n 
        norm_x = norm(R(j:end,j));
        s     = -sign(R(j,j));
        u_1    = R(j,j) - s*norm_x;
        w     = R(j:end,j)/u_1;
        w(1)  = 1;
        tau   = -s*u_1/norm_x;
        R(j:end,:) = R(j:end,:)-(tau*w)*(w'*R(j:end,:));
        Q(:,j:end) = Q(:,j:end)-(Q(:,j:end)*w)*(tau*w)';
        contador = contador + 3;
    end
    z = inv(Q)*b;
    sol = inv(R)*z; 
    error = norm(eye(m)-inv(Q*R)*A);
    contador = contador + 1;
end