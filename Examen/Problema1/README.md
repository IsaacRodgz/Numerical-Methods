# Problema1

## Uso

* Solucion sistema de ecuaciones:

$ make solve
$ ./runTest A_1000.mtx b_1000.vec > X_1000.vec

Se resuelve el sistema para matriz A_1000.mtx y vector b_1000.vec guardando la solucion en archivo X_1000.vec

* Calculo de inversa

$ make inverse
$ ./runTest A_1000.mtx > INV_1000.mtx

Se calcula la inversa de matriz A_1000.mtx guardando la matriz inversa en archivo INV_1000.mtx

* El archivo b_5000.vec se guarda con las dos dimensiones en el primer renglon debido a que la forma en que tengo implementados los algoritmos leen el vector como una matrix de nx1 como se había manejado en las tareas 

## Built With

* [C11] en Ubuntu 16.04 LTS

## Authors

* **Isaac Rodríguez Bribiesca**

