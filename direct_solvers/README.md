# Direct solvers for linear systems

Funciones para resolver sistemas de ecuaciones lineales (A*x = b) y determinantes:

* Diagonal
* Triangular Inferior
* Tringular superior
* Gauss sin pivote
* Gauss con pivoteo completo
* Factorización Doolittle
* Factorización Cholesky modificado

Así como la función para encontrar la matriz inversa.

## Estructura

Raiz:

* Biblioteca solve_matrix_direct.c, solve_matrix_direct.h con todas las funciones que resuelven los sistemas de ecuaciones, inversa y determinantes.

* Biblioteca matrix_struct.c, matrix_struct.h con funciones propias de una matriz como:

    - Leer matriz desde archivo
    - Imprimirla en pantalla
    - Multiplicar matrices
    - Verificar si 2 matrices son iguales, elemento a elemento
    
    Así como estructura Matriz con elementos:
        
    - Cols : Número de columnas de la matriz 
    - Rows : Número de renglones de la matriz
    - Data : Apuntador a double que guardará la matriz

* Carpetas con los nombres de los métodos. Cada carpeta incluye:

    - Archivo .c que ejecuta y prueba el método
    - Archivo makefile para compilar la prueba

## Usage

- Posicionarse dentro de cada carpeta "metodoSolver" y ejecutar comando "make" para compilar la prueba "metodo_solver_test.c"

- Se genera archivo ejeutable llamado "runTest"

- Ejecutar comando "./runTest /path/to/matrix.txt /path/to/vector.txt". Donde se pasan 2 argumentos, la ruta del archivo donde está la matriz y la ruta donde está el vector respectivamente. Excepto para el método de matriz inversa que sólo recibe un argumento correspondiente a la matriz.

## Example

Dentro carpeta "DoolittleSolver":

/path/DoolittleSolver >> make

/path/DoolittleSolver >> ./runTest ../MATRICES/M_SMALL.txt ../MATRICES/V_SMALL.txt

/path/DoolittleSolver >> make clean

## Test output

La ejecución de cada ejecutable runTest muestra en pantalla:

- Impresión en consola de las matrices leidas

- Vector solución x_solver al sistema A*x = b, calculada por el algoritmo

- Resultado de prueba de solución del algoritmo comparando los vectores b y A*x_solver. Dicha comparación se realiza con una tolerancia de 0.000000001

- Determinante de la matriz A leida

## Built With

* [C11] en Ubuntu 16.04 LTS

## Authors

* **Isaac Rodríguez Bribiesca**

