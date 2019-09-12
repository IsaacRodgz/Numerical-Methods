#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../matrix_struct.h"

void test_solve_jacobi(const char *matrix_filename){

    printf("\nRead matrix A");
    Matrix *A;
    A = read_matrix(matrix_filename, 1);
    //print_matrix(A);

    // Vector that will contain all the eigenvetors of A
    Matrix *eigenVects = malloc( sizeof( eigenVects ) );
    eigenVects->rows = A->rows;
    eigenVects->cols = A->cols;
    eigenVects->data = malloc( eigenVects->rows * eigenVects->cols * sizeof( eigenVects->data ) );

    int numIters = 10000;
    double epsilon = 0.000000001;;

    jacobiSolver(A, eigenVects, numIters, epsilon);

    print_matrix(A);
    print_matrix(eigenVects);
}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, /path/to/vector\n\n");
    }

    else
        test_solve_jacobi(argv[1]);

}
