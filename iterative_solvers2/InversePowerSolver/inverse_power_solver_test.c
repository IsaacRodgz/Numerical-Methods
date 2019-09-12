#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../matrix_struct.h"

void test_inverse_solve_power(const char *matrix_filename){

    printf("\nRead matrix A");
    Matrix *A;
    A = read_matrix(matrix_filename, 1);
    //print_matrix(A);

    // Vector that will contain the dominant eigenvetor corresponding to the dominant eigenvalue found
    Matrix *eigenVec = malloc( sizeof( eigenVec ) );
    eigenVec->rows = A->rows;
    eigenVec->cols = 1;
    eigenVec->data = malloc( eigenVec->rows * eigenVec->cols * sizeof( eigenVec->data ) );

    double eigenVal = 0;
    int numIters = 10000;
    double epsilon = 0.000000001;

    inversePowerSolver(A, eigenVec, &eigenVal, numIters, epsilon);

    printf("----------------------------------------------\n\n");
    printf("Dominant eigenvalue found: %lf\n", eigenVal);
    printf("\nDominant eigenvector found:\n");
    print_matrix( eigenVec );

    printf("----------------------------------------------\n\n");
    A = read_matrix(matrix_filename, 0);
    printf("Product of A*eigenvector = lambda*eigenvector: \n");
    print_matrix( multiply( A, eigenVec ) );

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, /path/to/vector\n\n");
    }

    else
        test_inverse_solve_power(argv[1]);

}
