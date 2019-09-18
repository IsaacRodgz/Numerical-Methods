#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative3.h"
#include "../matrix_struct.h"

void test_solve_cgradient(const char *matrix_filename, const char *vector_filename){

    Matrix *A;
    A = read_matrix(matrix_filename, 1);

    Matrix *b;
    b = read_matrix(vector_filename, 1);

    Matrix *x;
    x = read_matrix(vector_filename, 1);

    int numIters = 50;
    double epsilon = 0.000000001;

    cGradientSolver(A, b, x, numIters, epsilon);

    print_matrix(x);

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, /path/to/vector\n\n");
    }
    if(argc==2){

        printf("\nerror: One missing argument is required: /path/to/vector\n\n");
    } else
        test_solve_cgradient(argv[1], argv[2]);

}
