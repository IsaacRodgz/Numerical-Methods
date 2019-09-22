#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative3.h"
#include "../matrix_struct.h"

void test_solve_rayleigh(const char *matrix_filename){

    Matrix *A;
    A = read_matrix(matrix_filename, 1);

    Matrix *eigenVect = malloc( sizeof( eigenVect ) );
    eigenVect->rows = A->rows;
    eigenVect->cols = 1;
    eigenVect->data = malloc( eigenVect->rows * eigenVect->cols * sizeof( *eigenVect->data ) );

    double lambda = 0.000081;
    //double lam;
    //sscanf(lambdaInit, "%lf", &lambda);

    int numIters = 1000;
    double epsilon = 0.0000000001;

    printf("\nValor propio inicial: %e\n", lambda);
    rayleighSolver(A, eigenVect, &lambda, numIters, epsilon);
    printf("\nValor propio aproximado: %e\n\n", lambda);

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: One missing argument is required: /path/to/matrix initialLambda\n\n");
    }

    else{
        //double lambda = strtod(argv[2], NULL);
        test_solve_rayleigh(argv[1]);
    }

}
