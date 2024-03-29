#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative3.h"
#include "../matrix_struct.h"

void test_solve_qr(const char *matrix_filename){

    Matrix *A;
    A = read_matrix(matrix_filename, 1);

    // Matrix that will contain the k dominant eigenvetors corresponding to the dominant eigenvalues found
    Matrix *eigenVects = malloc( sizeof( eigenVects ) );
    eigenVects->rows = A->rows;
    eigenVects->cols = A->cols;
    eigenVects->data = malloc( eigenVects->rows * eigenVects->cols * sizeof( eigenVects->data ) );

    int numIters = 100;
    double epsilon = 0.000000001;

    QRSolve(A, eigenVects, numIters, epsilon);

    printf("\nEigenvalues:\n\n");
    for (int i = 0; i < A->rows; i++) {
        printf("%e\n", A->data[(A->cols+1)*i]);
    }

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: One missing argument is required: /path/to/matrix\n\n");
    } else

        test_solve_qr(argv[1]);

}
