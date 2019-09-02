#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../matrix_struct.h"
#include "../plot.h"

void test_solve_power(const char *matrix_filename){

    printf("\nRead matrix A");
    Matrix *A;
    A = read_matrix(matrix_filename, 1);
    print_matrix(A);

    // Vector that will contain the dominant eigenvetor corresponding to the dominant eigenvalue found
    Matrix *eigenVec = malloc( sizeof( eigenVec ) );
    eigenVec->rows = A->rows;
    eigenVec->cols = 1;
    eigenVec->data = malloc( eigenVec->rows * eigenVec->cols * sizeof( eigenVec->data ) );

    for (size_t i = 0; i < eigenVec->rows; i++) {
        eigenVec->data[i] = 1.0;
    }

    double eigenVal = 0.0;
    int numIters = 200;
    double epsilon = 0.0000001;

    double* lambdaTrace;

    lambdaTrace = powerSolver(A, eigenVec, &eigenVal, numIters, epsilon);

    printf("Dominant eigenvalue found: %lf\n", eigenVal);
    printf("\nDominant eigenvector found:\n");
    print_matrix( eigenVec );

    plotData( lambdaTrace, numIters );

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, /path/to/vector\n\n");
    }

    if(argc==2){

        test_solve_power(argv[1]);
    }

}
