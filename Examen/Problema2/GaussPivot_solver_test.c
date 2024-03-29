#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve_matrix_direct.h"
#include "matrix_struct.h"

void test_gaussPivot(const char *matrix_filename, const char *vector_filename){

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

    x_solve = solve_pivot(A, b);

    printf("-----------------------------------------------------------\n\n");

    //printf("\nSolution of system by Gauss with complete pivoting, x_solve:\n");
    //print_matrix(x_solve);

    A_temp = read_matrix(matrix_filename, 0);
    b = read_matrix(vector_filename, 0);

    printf("\nComparing A*x_solve with b...\n");

    double total = 0;

    for(int i = 0; i < A_temp->rows; i++){

        double sum = 0;

        for(int j = 0; j < A_temp->cols; j++){

            sum += A_temp->data[A_temp->cols*i + j]*x_solve->data[j];
        }

        total += (sum - b->data[i])*(sum - b->data[i]);
    }

    printf("\nError: %lf\n", total);
}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, /path/to/vector\n\n");
    }

    if(argc==2){

        printf("\nerror: One missing argument is required: /path/to/vector\n");
    }

    if(argc==3){

        test_gaussPivot(argv[1], argv[2]);
    }
}
