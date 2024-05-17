# Algoritmos Metodos Numericos




### Punto Fijo

```
function FIJO(x0, es, imax, gx)
    xr = x0;
    iter = 0;
    g = str2func(['@(x)' gx]);
    do = 0;
    fprintf(' iteracion Resultado error\n')
    while (do == 0)
        xrold = xr;
        xr = g(xrold);
        iter = iter + 1;
        if (xr ~= 0)
            ea = abs((xr - xrold) / xr) * 100;
        end
        fprintf(' %4.0f %4.4f %4.4f\n', iter, xr, ea)
        if ((ea < es) || (iter >= imax))
            break;
        end
    end
end

```

### Metodo de Bairstow
```

while true % Este ciclo nos permite validar que el primer coeficiente sea 1
    fx = input('Ingrese la matriz de coeficientes de la función fx: ');
    if fx(1) == 1
        break
    else
        disp('ERROR!!!!!')
        disp('El primer coeficiente debe ser igual a 1, por favor intente de nuevo.')
    end
end

% Solicitar al usuario ingresar los valores de r y s
r = input('Ingrese el valor de r: ');
s = input('Ingrese el valor de s: ');

format long; % Usaremos la mayor cantidad de dígitos posibles.
Tol = 0.0001;
n = length(fx); % Tamaño de la función
it = 0;
Maxit = 800; % Ponemos un máximo de iteraciónes en caso de que no converger
Roots = sparse(1, n - 1); % Creamos la matriz que nos almacenará las raices
cRoot = 0; % Contador de la matriz
Er = 100; % Para iniciar los errores los colocaremos en un 100%
Es = 100;

while true
    it = it + 1; % Se incrementa en 1 la iteración
    if (n - 1 >= 3)
        bn = fx(1);
        Mbi(1) = bn; % Mbi es una matriz que almacenará todos los bi
        bn1 = fx(2) + r * bn;
        Mbi(2) = bn1;
        i = 3;
        while i <= n
            if i == 3
                bi = fx(i) + r * bn1 + s * bn;
            else
                bi = fx(i) + r * Mbi(i - 1) + s * Mbi(i - 2);
            end
            Mbi(i) = bi;
            i = i + 1;
        end
        i = 3;
        cn = bn;
        Mci(1) = cn; % Mci es una matriz que almacenará todos los ci
        cn1 = bn1 + r * cn;
        Mci(2) = cn1;
        while i <= n - 1
            if i == 3
                ci = Mbi(3) + r * cn1 + s * cn;
            else
                ci = Mbi(i) + r * Mci(i - 1) + s * Mci(i - 2);
            end
            Mci(i) = ci;
            i = i + 1;
        end
        % Se resuelve el sistema de ec
        b0 = Mbi(n);
        b1 = Mbi(n - 1);
        c1 = Mci(n - 1);
        c2 = Mci(n - 2);
        c3 = Mci(n - 3);
        s1 = ((-b0 * c2 + c1 * b1) / (-c3 * c1 + c2 ^ 2));
        r1 = ((-b1 - c3 * s1) / c2);
        r = r + r1;
        s = s + s1;
        % Error
        Er = abs((r1 / r) * 100);
        Es = abs((s1 / s) * 100);
        if ((abs(Er)) && (abs(Es)) <= Tol) % Cuando llegamos al error tolerado se calculan las raices en esa iteración
            cRoot = cRoot + 1;
            x1 = (r + sqrt((r ^ 2) + 4 * s)) / 2;
            Roots(cRoot) = x1;
            x2 = (r - sqrt((r ^ 2) + 4 * s)) / 2;
            cRoot = cRoot + 1;
            Roots(cRoot) = x2;
            parar = true;
            c = length(Mbi);
            while parar
                if abs(Mbi(c)) >= 0.0001
                    break;
                else
                    c = c - 1;
                end
            end
            fx = sparse(1, 1);
            while c > 0
                fx(c) = Mbi(c); % Redefinimos fx
                c = c - 1;
            end
            n = length(fx);
        end
    else
        break
    end
end

if (n - 1) == 2 % Cuando la ecuación es de segundo grado se resuelve con el método de fórmula general
    it = it + 1; % Se incrementa en 1 la iteración
    disc = (Mbi(2) ^ 2) - 4 * bn * Mbi(3); % Se calcula el discriminante
    x1 = (-Mbi(2) + sqrt(disc)) / (2 * bn);
    cRoot = cRoot + 1;
    Roots(cRoot) = x1;
    x1 = (-Mbi(2) - sqrt(disc)) / (2 * bn);
    cRoot = cRoot + 1;
    Roots(cRoot) = x1;
elseif (n - 1) == 1 % Cuando la ecuación es de primer grado solamente se hace la sustitución.
    it = it + 1; % Se incrementa en 1 la iteración
    if (abs(fx(2)) >= Tol)
        cRoot = cRoot + 1;
        x1 = -fx(2) / fx(1);
    else
        x1 = 0;
    end
    Roots(cRoot) = x1;
    x1;
end

if (n - 1) == 0 % Si el valor de la ecuacion que se introduce es una constante
    disp('La ecuación no presenta incógnitas, es una constante igual a:')
    constante = fx;
end

Roots
Metodo de Bairstow - Matlab.txt
Mostrando Metodo de Bairstow - Matlab.txt.
```

### Metodo Newton-Rapson
```

syms x
fprintf('Metodo de Newton Rapson \n');
funcion = input('Ingrese la funcion f(x): ');   % recibe cualquier tipo de funcion
xi = input('Ingrese el valor de Xi: ');  % Guarda el Valor de X
n = input('Ingrese el valor de numero de Iteraciones: '); % Es el numero de Veces que se quiere iterar
tolerancia = input('Ingrese la tolerancia: '); % Solicita la tolerancia al usuario
fprintf('\n');
fprintf('-------------|---------------------|------------------------\n');
fprintf('Iteracion    |     Estimacion      |      Error aprox      |\n');
fprintf('-------------|---------------------|------------------------\n');
fxi = inline(funcion); % Guarda la funcion original para posterior ser evaluada con el valor de xi

z = diff(funcion,1); % deriva la funcion, por primera vez
d1 = inline(z); % Guarda la primera derivada de la funcion

a = diff(funcion,2); % deriva la funcion por segunda vez
d2 = inline(a); % Guarda la segunda derivada de la funcion

for i = 1: n
    funvalu = fxi(xi); % se Valua la funcion Original con el valor de Xi
    d1valuada = d1(xi); % evalua y guarda la funcion Derivada por primera vez con el valor que tiene xi
    d2valuada = d2(xi); % evalua y guarda la funcion Derivada por segunda vez con el valor que tiene xi
    Estimacion = xi-((funvalu*d1valuada)/((d1valuada)^2-funvalu*d2valuada)); % Calcula la Estimacion
    Error = abs((Estimacion - xi) / Estimacion)*100; % Calcula el Error de Aproximacion
    fprintf('    %2d       |    %12.9f     |    %12.8f       |\n', i, Estimacion, Error);
    fprintf('-------------|---------------------|------------------------\n');
    
    if tolerancia>=Error || i == n % Verifica si se alcanzó la tolerancia o el número máximo de iteraciones
        fprintf('Proceso Finalizado...\n');
        break; % Sale del bucle
    end
    xi = Estimacion; % Asigna el nuevo valor de xi para la siguiente iteración
end
Metodo_NewtonRaphson_Modificado.m
Mostrando Metodo_NewtonRaphson_Modificado.m.

```
### Metodo Falsa posicion 

```
function falsa_posicion_main()
    funcion_str = input('Ingrese la función f(x) en términos de x: ', 's');
    fx = str2func(['@(x) ', funcion_str]);
    p0 = input('Ingrese el valor inicial p0: ');
    p1 = input('Ingrese el valor inicial p1: ');
    TOL = input('Ingrese la tolerancia TOL: ');
    N0 = input('Ingrese el número máximo de iteraciones N0: ');
    
    if fx(p0) * fx(p1) >= 0
        error('Los valores iniciales p0 y p1 deben tener signos opuestos para aplicar la regla falsa.');
    end
    
    ta = zeros(N0, 7);
    i = 1;
    a = p0;
    b = p1;
    fa = fx(a);
    fb = fx(b);
    
    fprintf('Iteración   a        c        b        fa       fc       fb       Tramo\n');
    fprintf('-------------------------------------------------------------------------------\n');
    
    while i <= N0
        c = (a * fb - b * fa) / (fb - fa);
        fc = fx(c);
        ta(i, :) = [a, c, b, fa, fc, fb, b - a];
        
        fprintf('%d      %.6f   %.6f   %.6f   %.6f   %.6f   %.6f   %.6f\n', i, a, c, b, fa, fc, fb, b - a);
        
        if abs(fc) < TOL
            fprintf('\nRaíz aproximada: %.6f\n', c);
            fprintf('Error final: %.6f\n', b - a);
            return;
        end
        
        if fa * fc < 0
            b = c;
            fb = fc;
        else
            a = c;
            fa = fc;
        end
        
        i = i + 1;
    end
    
    fprintf('El método falló después de %d iteraciones.\n', N0);
end

```


### Metodo de la secante

```
% Script para resolver una función usando el método de la secante y mostrar una tabla de iteraciones con títulos

% Solicitar al usuario que ingrese la función
func_str = input('Introduce la función f(x): ', 's');
f = str2func(['@(x)', func_str]);

% Solicitar los valores iniciales
x0 = input('Introduce el valor inicial x0: ');
x1 = input('Introduce el valor inicial x1: ');

% Solicitar la tolerancia
tol = input('Introduce la tolerancia: ');

% Llamar a la función del método de la secante y obtener la tabla de iteraciones
[table, root] = secant_method(f, x0, x1, tol);

% Definir los títulos de las columnas
column_titles = {'Iteración', 'x_i-1', 'x_i', 'f(x_i-1)', 'f(x_i)', 'x_i+1', 'Error relativo'};

% Mostrar la tabla de iteraciones con títulos
disp('Tabla de iteraciones:');
disp(repmat('-', 1, 74)); % Línea horizontal
fprintf('%-10s%-10s%-10s%-12s%-12s%-12s%-15s\n', column_titles{:}); % Títulos de las columnas
disp(repmat('-', 1, 74)); % Línea horizontal
for i = 1:size(table, 1)
    fprintf('%-10d%-10.6f%-10.6f%-12.6f%-12.6f%-12.6f%-15.6e\n', table(i, :));
end
disp(repmat('-', 1, 74)); % Línea horizontal

% Mostrar el resultado final
fprintf('La raíz aproximada es: %f\n', root);

% Función del método de la secante
function [table, root] = secant_method(f, x0, x1, tol)
    % f: función a la cual queremos encontrar la raíz
    % x0, x1: dos puntos iniciales
    % tol: tolerancia para el criterio de convergencia
    
    % Inicializar la tabla de iteraciones
    table = zeros(1, 7);
    
    % Evaluar la función en los puntos iniciales
    f_x0 = f(x0);
    f_x1 = f(x1);
    
    % Inicializar el número de iteraciones
    iter = 0;
    
    % Iterar hasta alcanzar la tolerancia
    while abs(f_x1) > tol
        % Calcular la siguiente aproximación usando el método de la secante
        x2 = x1 - f_x1 * ((x1 - x0) / (f_x1 - f_x0));
        
        % Calcular el error relativo
        if x2 ~= 0
            relative_error = abs(x2 - x1) / abs(x2);
        else
            relative_error = abs(x2 - x1);
        end
        
        % Incrementar el contador de iteraciones
        iter = iter + 1;
        
        % Agregar la información de la iteración actual a la tabla
        table(iter, :) = [iter, x0, x1, f_x0, f_x1, x2, relative_error];
        
        % Actualizar los valores para la siguiente iteración
        x0 = x1;
        f_x0 = f_x1;
        x1 = x2;
        f_x1 = f(x1);
    end
    
    % La raíz aproximada es la última x1 calculada
    root = x1;
end

```


### Metodo Biseccion

```
function raiz = biseccion(f, a, b, tol)
    % Función que encuentra la raíz de una función utilizando el método de bisección
    
    fa = f(a);
    fb = f(b);
    
    if fa * fb > 0
        error('La función no cambia de signo en el intervalo dado.');
    end
    
    fprintf('Iteración |   a   |   b   |   c   |  f(c) \n');
    fprintf('-----------------------------------------\n');
    
    iter = 0;
    while (b - a) / 2 > tol
        c = (a + b) / 2;
        fc = f(c);
        
        fprintf('%5d     | %.4f | %.4f | %.4f | %.4f\n', iter, a, b, c, fc);
        
        if abs(fc) < tol
            raiz = c;
            return;
        end
        
        if fa * fc < 0
            b = c;
            fb = fc;
        else
            a = c;
            fa = fc;
        end
        
        iter = iter + 1;
    end
    
    raiz = (a + b) / 2;
end

% Ejemplo de uso:
f_string = input('Ingrese la función f(x): ', 's'); % Solicitamos la función al usuario como una cadena
f = str2func(['@(x)' f_string]); % Convertimos la cadena en una función anónima
a = input('Ingrese el valor de a: ');
b = input('Ingrese el valor de b: ');
tol = input('Ingrese la tolerancia deseada: ');

% Llamada a la función biseccion con los datos proporcionados
raiz = biseccion(f, a, b, tol);

fprintf('La raíz aproximada es: %.6f\n', raiz);

```


