#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../matrix_struct.h"

void test_solve_k_power(const char *matrix_filename, int numEigen){

    printf("\nRead matrix A");
    Matrix *A;
    A = read_matrix(matrix_filename, 1);

    // parameter to define number of eigenvectors and eigenvalues to find
    int k = numEigen;

    // Matrix that will contain the k dominant eigenvetors corresponding to the dominant eigenvalues found
    Matrix *eigenVects = malloc( sizeof( eigenVects ) );
    eigenVects->rows = k;
    eigenVects->cols = A->rows;
    eigenVects->data = malloc( eigenVects->rows * eigenVects->cols * sizeof( eigenVects->data ) );

    // Vector that will contain the k dominant eigenvalues corresponding to the dominant eigenvectors found
    Matrix *eigenVals = malloc( sizeof( eigenVals ) );
    eigenVals->rows = k;
    eigenVals->cols = 1;
    eigenVals->data = malloc( eigenVals->rows * eigenVals->cols * sizeof( eigenVals->data ) );

    // Parameters for iteration
    int numIters = 10000;
    double epsilon = 0.000000001;

    // Find eigenvetors and eigenvalues

    kPowerSolver(A, eigenVects, eigenVals, numIters, epsilon, k);

    printf("----------------------------------------------\n\n");
    printf("Dominant eigenvalues found:\n");
    print_matrix( eigenVals );
    printf("\nDominant eigenvectors found:\n");
    print_matrix( eigenVects );

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, numEigen\n\n");
    }
    if(argc==2){

        printf("\nerror: One missing argument is required: numEigen\n\n");
    }

    else
        test_solve_k_power(argv[1], atoi(argv[2]));

}
