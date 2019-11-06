#include <stdio.h>
#include <stdlib.h>
#include "solve_matrix_direct.h"
#include "matrix_struct.h"

void test_solveInverse(const char *matrix_filename){

    Matrix *A;
    Matrix *A_solve;

    //printf("\nRead matrix A");
    A = read_matrix(matrix_filename, 1);
    //print_matrix(A);

    A_solve = solve_inverse(A);

    print_matrix(A_solve);

    /*
    printf("-----------------------------------------------------------\n\n");

    // Read A again because solve_no_pivot() overwrites it

    A = read_matrix(matrix_filename, 0);

    // Build Identity matrix I

    Matrix *I = malloc( sizeof( I ) );

    I->rows = A->rows;
    I->cols = A->cols;
    I->data = malloc( I->rows*I->cols*sizeof( I->data ) );

    for(int i = 0; i < I->rows; i++){

        for(int j = 0; j < I->cols; j++){

            if( i == j )
                I->data[ i*I->cols + j ] = 1.0;
            else
                I->data[ i*I->cols + j ] = 0.0;
        }
    }


    printf("\nComparing A * A^-1 with identity matrix I...\n");
    if(equals(I, multiply( A_solve, A ), 0.000000001) == 0)
        printf("\nTest FAILED. A * A^-1 != I for given A and b.\n\n");
    else
        printf("\nTest PASSED. A * A^-1 = I for given A and b.\n\n");
    */
}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: One missing argument is required: /path/to/vector\n");
    }

    if(argc==2){

        test_solveInverse(argv[1]);
    }
}
