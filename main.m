clear all

A = load('sistemas/289x289/A289.dat');
b = load('sistemas/289x289/b289.dat');

%A = load('sistemas/1089x1089/A1089.dat');
%b = load('sistemas/1089x1089/b1089.dat');

%A = load('sistemas/4225x4225/A4225.dat');
%b = load('sistemas/4225x4225/b4225.dat');

%[valuesToGraphX,valuesToGraphF,valuesToGraphError, iteraciones] = multivariable_newton();

%[sola, errora] = householder(A, b);
%[solb, errorb] = gauss_jacobi(A, b);
%[solc, errorc] = gauss_seidel(A, b);

[n,m] = size(A);



%Iterativos para matrices esparcidas y grandes 
%Redondea, no trucan
%if(SparseMatrix(A) && n==m && n > 1000)
%    if(DominantDiagonal(A))
        [sol_gs, error_gs] = gauss_seidel(A,b);
%    end
    
%    if(DominantDiagonal(A))
        [sol_gj, error_gj] = gauss_jacobi(A,b);
%    end   

%    if(PositiveMatrix(A) &&  SimetricMatrix(A)) %Si la matriz es definida positiva y simétrica, se puede usar cholesky O(n^3/3)
        [sol_chol, error_chol] = cholesky(A,b);
%    end
%    if(det(A) ~= 0) %Si tiene inversa se puede usar doolite O(2n^3/3)
        [sol_dool, error_dool] = doolittle(A,b);
%    end

%    
%end

%Este método se suele usar cuando la estabilidad numérica es prioritaria
%(se quiere minimizar el error) O(3n^2 (m - n/3))
%if(SparseMatrix(A) && (m>=n))
    [sol_giv, error_giv] = givens(A,b);
%end

%if(~SparseMatrix(A) && m>=n)
    [sol_hous, error_hous] = householder(A,b); 
%end

%O(2*m*n^2)
%if(m>n)
    [sol_QR, error_QR] = QR(A,b);
%end



%figure
%plot(error_gs)

%figure
%plot(error_gj)

%SimetricMatrix(matrix);

%DominantDiagonal(matrix);

%PositiveMatrix(A);

%SemiPositiveMatrix(matrix2);

%NegativeMatrix(matrix3);