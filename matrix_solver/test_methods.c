#include <stdio.h>
#include <stdlib.h>
#include "solve_matrix_direct.h"
#include "matrix_struct.h"


void test_read_matrix(const char *matrix_filename, const char *vector_filename){

    Matrix *matrix;
    Matrix *vector;

    matrix = read_matrix(matrix_filename, 1);
    print_matrix(matrix);

    vector = read_matrix(vector_filename, 1);
    print_matrix(vector);
}

void test_mult(const char *matrix_filename, const char *vector_filename){

    Matrix *x = malloc( sizeof( x ) );
    Matrix *y = malloc( sizeof( y ) );
    Matrix *z = malloc( sizeof( z ) );

    x->rows = 3;
    x->cols = 3;
    x->data = malloc( x->rows*x->cols*sizeof( x->data ) );
    double data1[] = {1,2,3,4,5,6,1,1,1};
    x->data = data1;

    print_matrix(x);

    y->rows = 4;
    y->cols = 1;
    y->data = malloc( y->rows*y->cols*sizeof( y->data ) );
    double data2[] = {1,0,1,1};
    y->data = data2;

    print_matrix(y);

    z = multiply(x, y);

    print_matrix(z);
}

void test_diagonal_solve(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    Matrix *b;
    Matrix *x_solve;

    A = read_matrix(matrix_filename, 1);
    print_matrix(A);

    b = read_matrix(vector_filename, 1);
    print_matrix(b);

    x_solve = solve_diagonal(A, b);

    printf("\nSolution of diagonal system, x_solve:\n");
    print_matrix(x_solve);

    if(equals(b, multiply(A, x_solve), 0.00000000001) == 0)
        printf("\nTest FAILED. A * x_solve != b for given A and b.\n");
    else
        printf("\nTest PASSED. A * x_solve = b for given A and b.\n");
}

void test_lower_triag(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    Matrix *b;
    Matrix *x_solve;

    A = read_matrix(matrix_filename, 1);
    print_matrix(A);

    b = read_matrix(vector_filename, 1);
    print_matrix(b);

    x_solve = solve_lower_triang(A, b, 0);

    printf("\nSolution of lower triangular system, x_solve:\n");
    print_matrix(x_solve);

    if(equals(b, multiply(A, x_solve), 0.00000000001) == 0)
        printf("\nTest FAILED. A * x_solve != b for given A and b.\n");
    else
        printf("\nTest PASSED. A * x_solve = b for given A and b.\n");
}

void test_upper_triag(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    Matrix *b;
    Matrix *x_solve;

    A = read_matrix(matrix_filename, 1);
    print_matrix(A);

    b = read_matrix(vector_filename, 1);
    print_matrix(b);

    x_solve = solve_upper_triang(A, b, 0);

    printf("\nSolution of lower triangular system, x_solve:\n");
    print_matrix(x_solve);

    if(equals(b, multiply(A, x_solve), 0.00000000001) == 0)
        printf("\nTest FAILED. A * x_solve != b for given A and b.\n");
    else
        printf("\nTest PASSED. A * x_solve = b for given A and b.\n");
}

void test_solve_no_pivot(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    Matrix *b;
    Matrix *x_solve;

    A = read_matrix(matrix_filename, 1);
    print_matrix(A);

    b = read_matrix(vector_filename, 1);
    print_matrix(b);

    x_solve = solve_no_pivot(A, b);

    printf("\nSolution vector x:\n");
    print_matrix(x_solve);

    A = read_matrix(matrix_filename, 0);
    b = read_matrix(vector_filename, 0);

    if(equals(b, multiply(A, x_solve), 0.00000000001) == 0)
        printf("\nTest FAILED. A * x_solve != b for given A and b.\n");
    else
        printf("\nTest PASSED. A * x_solve = b for given A and b.\n");
}

void test_solve_pivot(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    Matrix *b;
    Matrix *x_solve;

    A = read_matrix(matrix_filename, 1);
    print_matrix(A);

    b = read_matrix(vector_filename, 1);
    print_matrix(b);

    x_solve = solve_pivot(A, b);

    printf("\nSolution vector x:\n");
    print_matrix(x_solve);

    A = read_matrix(matrix_filename, 0);
    b = read_matrix(vector_filename, 0);

    if(equals(b, multiply(A, x_solve), 0.00000000001) == 0)
        printf("\nTest FAILED. A * x_solve != b for given A and b.\n");
    else
        printf("\nTest PASSED. A * x_solve = b for given A and b.\n");
}

void test_solve_doolittle(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    Matrix *b;
    Matrix *x_solve;

    A = read_matrix(matrix_filename, 1);
    print_matrix(A);

    b = read_matrix(vector_filename, 1);
    print_matrix(b);

    x_solve = solve_doolittle(A, b, 1);

    printf("\nSolution vector x:\n");
    print_matrix(x_solve);

    A = read_matrix(matrix_filename, 0);
    b = read_matrix(vector_filename, 0);

    if(equals(b, multiply(A, x_solve), 0.00000000001) == 0)
        printf("\nTest FAILED. A * x_solve != b for given A and b.\n");
    else
        printf("\nTest PASSED. A * x_solve = b for given A and b.\n");

}

void test_solve_inverse(const char *matrix_filename){

    Matrix *A;
    Matrix *A_solve;

    A = read_matrix(matrix_filename, 1);
    print_matrix(A);

    A_solve = solve_inverse(A);

    A = read_matrix(matrix_filename, 0);

    print_matrix( multiply( A_solve, A ) );

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

    if(equals(I, multiply( A_solve, A ), 0.00000000001) == 0)
        printf("\nTest FAILED. A * x_solve != b for given A and b.\n");
    else
        printf("\nTest PASSED. A * x_solve = b for given A and b.\n");

}

void test_solve_cholesky_modified(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    Matrix *b;
    Matrix *x_solve;

    A = read_matrix(matrix_filename, 1);
    print_matrix(A);

    b = read_matrix(vector_filename, 1);
    print_matrix(b);

    x_solve = solve_cholesky_modified(A, b);

    //print_matrix(A);

    printf("\nSolution vector x:\n");
    print_matrix(x_solve);

    A = read_matrix(matrix_filename, 0);
    b = read_matrix(vector_filename, 0);

    if(equals(b, multiply(A, x_solve), 0.00000000001) == 0)
        printf("\nTest FAILED. A * x_solve != b for given A and b.\n");
    else
        printf("\nTest PASSED. A * x_solve = b for given A and b.\n");

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nUsage: gcc solve_matrix_direct.c test_methods.c matrix_path vector_path");
        printf("\nerror: the following arguments are required: matrix_path, vector_path\n");
    }

    if(argc==2){

        //printf("\nUsage: gcc solve_matrix_direct.c test_methods.c matrix_path vector_path");
        //printf("\nerror: the following arguments are required: vector_path\n");
        test_solve_inverse( argv[1] );
    }

    if(argc==3){

        test_solve_cholesky_modified(argv[1], argv[2]);
    }
}
