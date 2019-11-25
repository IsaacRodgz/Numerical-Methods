# Least squares and FEM

## Usage

- Posicionarse dentro de cada carpeta de los programas y ejecutar el siguiente comando:
    * make && ./test

- Para modelo lineal se usa archivo "lineal_data.txt"
- Para modelo cuadrático se usa archivo "quadratic_data.txt"

## Output

* Para modelo lineal y cuadrático se muestra por consola el modelo ajustado
* Para el método de elemento finito se muestran por consola parejas de datos (x_i, Phi_i),
  uno por cada elemento. Donde x_i representa el punto x que se usó como nodo y Phi_i representa
  el coeficiente ajustado para el elemento i.

## Built With

* [C11] en Ubuntu 16.04 LTS

## Authors

* **Isaac Rodríguez Bribiesca**
