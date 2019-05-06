function [sol,error,contador] = givens(A,b)

    % sol es vector solución
    % error es el error final
    % A matriz del sistema
    % b vector del sistema
    
    contador = 0;
    [m,n] = size(A);
    Q = eye(m);
    R = A;
    
    for i=1:n
        for k=i+1:m
            if (R(k,i) ~= 0)
                contador = contador + 1;
                
                raiz = sqrt(R(k,i)^2 + R(i,i)^2);
                s = -R(k,i)/raiz;
                c = R(i,i)/raiz;
                G = eye(m); 
                
                G(k,k) = c;
                G(i,i) = c;
                
                G(i,k) = s;
                G(k,i) = -s;
                
                Q = Q*G; 
                R = G'*R; 
                
            end
        end
        
        contador = contador + 1;
    end
    
      Y = inv(Q)*b;  
      sol = inv(R)*Y;
      error = norm(eye(n)-inv(Q*R)*A);
      contador = contador + 1;
end
