# Eigensolvers

* Rayleigh
* Iteración en subespacio
* QR

# Solver

* Gradiente conjugado

## Usage

- Posicionarse dentro de cada carpeta "metodoSolver" y ejecutar comando "make" para compilar el programa

- Se genera archivo ejeutable llamado "runTest"

- Ejecutar comando "./runTest rutaMatriz

## Arguments

* Rayleigh: rutaMatriz
* Subespacio: rutaMatriz numEigenVals
* QR: rutaMatriz
* Gradiente conjugado: rutaMatriz rutaVector

## Example

Dentro carpeta "SubspaceSolver":

/path/SubspaceSolver >> make

/path/SubspaceSolver >> ./runTest ../M_BIG.txt 20

/path/SubspaceSolver >> make clean

## Built With

* [C11] en Ubuntu 16.04 LTS

## Authors

* **Isaac Rodríguez Bribiesca**
