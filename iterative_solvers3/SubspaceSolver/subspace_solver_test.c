#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative3.h"
#include "../matrix_struct.h"

void test_solve_subspace(const char *matrix_filename, int numEigen){

    Matrix *A;
    A = read_matrix(matrix_filename, 1);

    // m eigenvectors
    // A[n][n]*FI[n][m] = FI[n][m]*L[m][m] -> Power method
    // FI[m][n]*A[n][n]*FI[n][m] = B[m][m] -> multiply() and test isDiagonal(B)
    // B[m][m] = Q[m][m]*L[m][m]*Q[m][m] -> Jacobi
    // FI[n][m]*Q[m][m] = FI[n][m] -> multiply()

    // parameter to define number of eigenvectors and eigenvalues to find
    int k = numEigen;

    // Matrix that will contain the k dominant eigenvetors corresponding to the dominant eigenvalues found
    Matrix *eigenVects = malloc( sizeof( eigenVects ) );
    eigenVects->rows = A->rows;
    eigenVects->cols = k;
    eigenVects->data = malloc( eigenVects->rows * eigenVects->cols * sizeof( eigenVects->data ) );

    // Vector that will contain the k dominant eigenvalues corresponding to the dominant eigenvectors found
    Matrix *eigenVals = malloc( sizeof( eigenVals ) );
    eigenVals->rows = k;
    eigenVals->cols = k;
    eigenVals->data = malloc( eigenVals->rows * eigenVals->cols * sizeof( eigenVals->data ) );

    int numIters = 200;
    double epsilon = 0.00000001;

    subspaceSolver(A, eigenVects, eigenVals, numIters, epsilon, k);

    printf("----------------------------------------------\n\n");
    printf("Eeigenvalues found:\n");
    for (int i = 0; i < eigenVals->rows; i++) {
        printf("%e\n", eigenVals->data[(eigenVals->cols+1)*i]);
    }
}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, numEigen\n\n");
    }
    if(argc==2){

        printf("\nerror: One missing argument is required: numEigen\n\n");
    } else
        test_solve_subspace(argv[1], atoi(argv[2]));

}
