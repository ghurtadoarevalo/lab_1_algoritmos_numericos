
format long
clear all
%close all 

%A = load('sistemas/289x289/A289.dat');
%b = load('sistemas/289x289/b289.dat');

%A = load('sistemas/1089x1089/A1089.dat');
%b = load('sistemas/1089x1089/b1089.dat');
    
A = load('sistemas/4225x4225/A4225.dat');
b = load('sistemas/4225x4225/b4225.dat');

estabilidadNumerica = false;
gestor = true;
graficar = false;
calcularTiempoMetodos = true;
[n,m] = size(A);

[valuesToGraphX,valuesToGraphF,valuesToGraphError, iteraciones] = multivariable_newton();

if(gestor == true)
    tic
    if(estabilidadNumerica == false)   
        %Iterativos para matrices esparcidas y grandes 
        %Redondea, no trucan
        if(SparseMatrix(A) && n==m && n >= 2000)
            if(DominantDiagonal(A) || (SimetricMatrix(A) && PositiveMatrix(A)))
                [sol_gs, error_gs, costo_gs] = gauss_seidel(A,b); 
            else
                [sol_gj, error_gj, costo_gj] = gauss_jacobi(A,b);
            end     
        elseif(n < 2000)
             %O(2*m*n^2)
             if(m>n)
                [sol_QR, error_QR, costo_QR] = QR(A,b);
             end
             if(n == m)     
                 if(PositiveMatrix(A) &&  SimetricMatrix(A)) %Si la matriz es definida positiva y simétrica, se puede usar cholesky O(n^3/3)
                    [sol_chol, error_chol, costo_chol] = cholesky(A,b);
                else           %O(2n^3/3)
                    [sol_dool, error_dool, costo_dool] = doolittle(A,b);
                end
             end
        end

    %Este método se suele usar cuando la estabilidad numérica es prioritaria
    %(se quiere minimizar el error) O(3n^2 (m - n/3))    
    else
        if(SparseMatrix(A))
            [sol_giv, error_giv, costo_giv] = givens(A,b);
        else(~SparseMatrix(A))
            [sol_hous, error_hous, costo_hous] = householder(A,b); 
        end  
    end
    tiempo_gestor = toc;

elseif(gestor == false && calcularTiempoMetodos == true)
    
    tic
    [sol_gs, error_gs, costo_gs] = gauss_seidel(A,b);
    tiempo_gs = toc;
    tic
    [sol_gj, error_gj, costo_gj] = gauss_jacobi(A,b);
    tiempo_gj = toc;
    tic
    [sol_QR, error_QR, costo_QR] = QR(A,b);
    tiempo_QR = toc;
    tic
    [sol_chol, error_chol, costo_chol] = cholesky(A,b);
    tiempo_chol = toc;
    tic
    [sol_dool, error_dool, costo_dool] = doolittle(A,b);
    tiempo_dool = toc;
    tic
    %[sol_giv, error_giv, costo_giv] = givens(A,b);
    error_giv = 0;
    costo_giv = 0;
    tiempo_giv = toc;
    tic
    [sol_hous, error_hous, costo_hous] = householder(A,b);
    tiempo_hous = toc;
    
    if(graficar == true)
        
        tiempo_graph = [];
        tiempo_graph = [tiempo_graph, tiempo_chol ];
        tiempo_graph = [tiempo_graph, tiempo_dool ];
        tiempo_graph = [tiempo_graph, tiempo_giv];
        tiempo_graph = [tiempo_graph, tiempo_gj];
        tiempo_graph = [tiempo_graph, tiempo_gs];
        tiempo_graph = [tiempo_graph, tiempo_hous];
        tiempo_graph = [tiempo_graph, tiempo_QR];
        
        plotTiempos(tiempo_graph)

        
        costo_graph = [];
        costo_graph = [costo_graph, costo_chol ];
        costo_graph = [costo_graph, costo_dool ];
        costo_graph = [costo_graph, costo_giv];
        costo_graph = [costo_graph, costo_gj];
        costo_graph = [costo_graph, costo_gs];
        costo_graph = [costo_graph, costo_hous];
        costo_graph = [costo_graph, costo_QR];

        plotCostos(costo_graph)

        err_graph = [];
        err_graph = [err_graph, error_chol];
        err_graph = [err_graph, error_dool ];
        err_graph = [err_graph, error_giv];
        err_graph = [err_graph, error_gj(end)];
        err_graph = [err_graph, error_gs(end)];
        err_graph = [err_graph, error_hous];
        err_graph = [err_graph, error_QR];

        %plot para todos los errores de los métodos
        plotErrors(err_graph)
    
    end

elseif(gestor == false && calcularTiempoMetodos == false)
    tic
    [sol_gs, error_gs, costo_gs] = gauss_seidel(A,b);
    [sol_gj, error_gj, costo_gj] = gauss_jacobi(A,b);
    [sol_QR, error_QR, costo_QR] = QR(A,b);
    [sol_chol, error_chol, costo_chol] = cholesky(A,b);
    [sol_dool, error_dool, costo_dool] = doolittle(A,b);
    [sol_giv, error_giv, costo_giv] = givens(A,b);
    [sol_hous, error_hous, costo_hous] = householder(A,b);
    tiempo_sin_gestor = toc;
    
    if(graficar == true)
        
        tiempo_graph = [];
        tiempo_graph = [tiempo_graph, tiempo_chol ];
        tiempo_graph = [tiempo_graph, tiempo_dool ];
        tiempo_graph = [tiempo_graph, tiempo_giv];
        tiempo_graph = [tiempo_graph, tiempo_gj];
        tiempo_graph = [tiempo_graph, tiempo_gs];
        tiempo_graph = [tiempo_graph, tiempo_hous];
        tiempo_graph = [tiempo_graph, tiempo_QR];
        
        plotTiempos(tiempo_graph)

        
        costo_graph = [];
        costo_graph = [costo_graph, costo_chol ];
        costo_graph = [costo_graph, costo_dool ];
        costo_graph = [costo_graph, costo_giv];
        costo_graph = [costo_graph, costo_gj];
        costo_graph = [costo_graph, costo_gs];
        costo_graph = [costo_graph, costo_hous];
        costo_graph = [costo_graph, costo_QR];

        plotCostos(costo_graph)

        sol_graph = [];
        sol_graph = [sol_graph, sol_chol ];
        sol_graph = [sol_graph, sol_dool ];
        sol_graph = [sol_graph, sol_giv];
        sol_graph = [sol_graph, sol_gj];
        sol_graph = [sol_graph, sol_gs];
        sol_graph = [sol_graph, sol_hous];
        sol_graph = [sol_graph, sol_QR];

        plotSolutions(sol_graph)

        err_graph = [];
        err_graph = [err_graph, error_chol];
        err_graph = [err_graph, error_dool ];
        err_graph = [err_graph, error_giv];
        err_graph = [err_graph, error_gj(end)];
        err_graph = [err_graph, error_gs(end)];
        err_graph = [err_graph, error_hous];
        err_graph = [err_graph, error_QR];

        %plot para todos los errores de los métodos
        plotErrors(err_graph)
    end
end



%plot para el error de seidel y jacobi
%plotErrorSeidelJacobi(error_gs,error_gj)

%SimetricMatrix(A);

%DominantDiagonal(matrix);

%PositiveMatrix(A);

%SemiPositiveMatrix(matrix2);

%NegativeMatrix(matrix3);