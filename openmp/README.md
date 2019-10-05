# Paralelización

## Usage

- Posicionarse dentro de cada carpeta de los programas y compilar con con la opción -fopenmp

- Ejecutar programa con argumento de número de hilos "./runTest numHilos

## Arguments

* Número de hilos a usar (En una máquina de k cores, si el argumento está entre 1 y k, el parámetro es equivalente al número de cores a usar)

## Example

Suma de vectores usando 4 hilos

Dentro carpeta "suma_vectores":

/path/suma_vectores >> gcc -fopenmp sum_vectors.c -o runTest

/path/suma_vectores >> ./runTest 4

## Built With

* [C11] en Ubuntu 16.04 LTS

## Authors

* **Isaac Rodríguez Bribiesca**
