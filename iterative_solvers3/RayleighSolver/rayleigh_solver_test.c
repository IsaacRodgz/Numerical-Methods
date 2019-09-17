#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative3.h"
#include "../matrix_struct.h"

void test_solve_rayleigh(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    A = read_matrix(matrix_filename, 1);

    Matrix *eigenVect;
    eigenVect = read_matrix(vector_filename, 1);

    int numIters = 1000;
    double epsilon = 0.00001;
    double lambda = -6.892;

    rayleighSolver(A, eigenVect, &lambda, numIters, epsilon);

    printf("\nValor propio aproximado: %e\n\n", lambda);

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix /path/to/vector\n\n");
    }
    else if (argc==2) {
        printf("\nerror: One missing argument is required: /path/to/vector\n\n");
    }

    else
        test_solve_rayleigh(argv[1], argv[2]);

}
