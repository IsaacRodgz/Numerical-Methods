#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../matrix_struct.h"

void test_solve_power(const char *matrix_filename){

    Matrix *A;
    printf("\nRead matrix A");
    A = read_matrix(matrix_filename, 1);
    print_matrix(A);


    Matrix *eigenVec = malloc( sizeof( eigenVec ) );
    eigenVec->rows = A->rows;
    eigenVec->cols = 1;
    eigenVec->data = malloc( eigenVec->rows * eigenVec->cols * sizeof( eigenVec->data ) );

    double eigenVal = 0.0;

    powerSolver(A, eigenVec, &eigenVal, 1000, 0.000000001);

    printf("Eigenvalue found: %lf\n\n", eigenVal);

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, /path/to/vector\n\n");
    }

    if(argc==2){

        test_solve_power(argv[1]);
    }

}
