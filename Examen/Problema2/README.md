# Problema 2

## Uso

* Gradiente Conjugado con matriz rala CSR

$ g++ -std=c++11 -O2 -fopenmp conjGradient.cpp solve_conjGradient.cpp -o runTest -O2 && ./runTest A_5000.mtx b_5000.vec

Para guardar solución x en archivo:
$ g++ -std=c++11 -O2 -fopenmp conjGradient.cpp solve_conjGradient.cpp -o runTest -O2 && ./runTest A_5000.mtx b_5000.vec > x_solve

* Gauss, Doolittle y Cholesky

$ make gauss
$ ./runTest A_5000.txt b_5000.vec

$ make doolittle
$ ./runTest A_5000.txt b_5000.vec

$ make cholesky
$ ./runTest A_5000.txt b_5000.vec

* El archivo b_5000.vec se guarda con las dos dimensiones en el primer renglon debido a que la forma en que tengo implementados los algoritmos leen el vector como una matrix de nx1 como se había manejado en las tareas 

## Built With

* [C11] en Ubuntu 16.04 LTS

## Authors

* **Isaac Rodríguez Bribiesca**

