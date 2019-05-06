function [valuesToGraphX,valuesToGraphF,valuesToGraphError, iteraciones] = multivariable_newton()
    syms x y z
    functions = [x^2+y-37;x-y^2-5;x+y+z-3];
    x0=[0; 0; 0];
    jac = jacobian(functions);
    inline_functions=inline(functions);
    inline_jac=inline(jac);
    error = norm(inline_functions(x0(1),x0(2),x0(3)),2);
    tolerancia=0.0000005;
    iteraciones = 30;
    i = 1;

    valuesToGraphX = [];
    valuesToGraphF = [];
    valuesToGraphError = [];
    valuesToGraphF = [valuesToGraphF,inline_functions(x0(1),x0(2),x0(3))];
    valuesToGraphX = [valuesToGraphX,x0];
    valuesToGraphError = [valuesToGraphError,error];

    while error >= tolerancia

        function_x = inline_functions(x0(1),x0(2),x0(3));
        jac_x = inline_jac(x0(1),x0(2));
        x1 = x0 - inv(jac_x) * function_x;
        fx1= inline_functions(x1(1),x1(2),x1(3));
        error = norm((fx1),2);

        valuesToGraphError = [valuesToGraphError,error];
        valuesToGraphX = [valuesToGraphX, x1];
        valuesToGraphF = [valuesToGraphF, fx1];

        if i > iteraciones
            fprintf(' No hay más iteraciones \n');
            break;
        end
        x0=x1;
        i=i+1;
    end
end

%figure
%plot(valuesToGraphX')
%figure
%plot(valuesToGraphF')
%figure
%plot(valuesToGraphError')

