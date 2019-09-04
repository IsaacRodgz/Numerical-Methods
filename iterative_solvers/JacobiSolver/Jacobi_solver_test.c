#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../matrix_struct.h"

void test_solve_jacobi(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    Matrix *A_temp;
    Matrix *b;
    Matrix *x_solve;

    printf("\nRead matrix A");
    A = read_matrix(matrix_filename, 1);
    //print_matrix(A);

    printf("\nRead vector b");
    b = read_matrix(vector_filename, 1);
    //print_matrix(b);

    x_solve = jacobiSolver(A, b, 5000, 0.000000001);

    //printf("\nSolution of system by Jacobi, x_solve:\n");
    //print_matrix(x_solve);

    printf("\nComparing A*x_solve with b...\n");
    //print_matrix(multiply(A_temp, x_solve));
    if(equals(b, multiply(A, x_solve), 0.0000000001) == 0)
        printf("\nTest FAILED. A * x_solve != b for given A and b.\n\n");
    else
        printf("\nTest PASSED. A * x_solve = b for given A and b.\n\n");

    printf("-----------------------------------------------------------\n\n");
}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, /path/to/vector\n\n");
    }

    if(argc==2){

        printf("\nerror: One missing argument is required: /path/to/vector\n");
    }

    if(argc==3){

        test_solve_jacobi(argv[1], argv[2]);
    }
}
