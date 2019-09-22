
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_struct.h"
#define TRUE 1
#define FALSE 0

Matrix * read_matrix(const char *filename, int verbose){

    FILE *file = fopen( filename, "r");

    if(file == NULL){

        perror(filename);

        return NULL;
    }

    // rows of Matrix
    int n;

    // columns of Matrix
    int m;

    // Read first line of file which contains num. rows and num. colums
    fscanf(file, "%d %d", &n, &m);

    // Matrix to fill with file
    Matrix *A = malloc( sizeof( A ) );

    // Initialize Matrix struct
    A->rows = n;
    A->cols = m;
    A->data = malloc( n*m*sizeof( A->data ) );

    int i,j;

    if(A->cols == 1){

        for(i = 0; i < n; i++){

            fscanf(file, "%lf ", (A->data + i));
        }
    }
    else{

        for(i = 0; i < n; i++){

            for(j = 0; j < m-1; j++){

                fscanf(file, "%lf ", (A->data + m*i + j));
            }

            fscanf(file, "%lf", (A->data + m*i + j));
        }
    }

    if(verbose != 0)
        printf("\nSuccesfully read: %s\n", filename);

    return A;
}

void print_matrix(Matrix *A){

    printf("\n");

    int rows = A->rows;
    int cols = A->cols;

    for(int i = 0; i < rows; i++){

        for(int j = 0; j < cols; j++){

            printf(" %lf", *(A->data + cols*i + j));
        }
        printf("\n");
    }
    printf("\n");
}

Matrix * multiply(Matrix *x, Matrix *y){

    if(x->cols != y->rows){

        printf("Error. Matrices dimensions don't match for multiplication\n");
        exit(-1);
        return NULL;
    }

    // z = x*y
    Matrix *z = malloc( sizeof( z ) );

    // Initialize Matrix struct
    z->rows = x->rows;
    z->cols = y->cols;
    z->data = malloc( z->rows*z->cols*sizeof( z->data ) );

    double sum;

    for(int i = 0; i < z->rows; i++){

        for(int j = 0; j < z->cols; j++){

            sum = 0;

            for(int k = 0; k < y->rows; k++){

                sum += x->data[ x->cols*i + k ] * y->data[ y->cols*k + j ];
            }

            z->data[ z->cols*i + j ] = sum;
        }
    }

    return z;
}

// Compare element-wise Matrix x and x_solve with tolerance epsilon

int equals(Matrix *x, Matrix *x_solve, double epsilon){

    int rows = x->rows;
    int cols = x->cols;

    double sum = 0;
    double diff;

    for(int i = 0; i < rows; i++){

        diff = ( x->data[i] - x_solve->data[i] );
        sum += diff * diff;
    }

    if( sqrt( sum ) < epsilon )
        return 1;
    else
        return 0;
}

// Copy content of vector x into y

void copy(Matrix *x, Matrix *y){

    int rows = x->rows;
    int cols = x->cols;

    for(int i = 0; i < rows; i++){

        for(int j = 0; j < cols; j++){

            y->data[ i*cols + j ] = x->data[ i*cols + j ];
        }
    }
}

// Helper function to verify if A is symmetric

int is_simetric(Matrix *A){

    for(int i = 0; i < A->rows; i++){

        for(int j = 0; j < A->cols; j++){

            if(A->data[ i*A->rows + j ] != A->data[ j*A->rows + i ])
                return FALSE;
        }
    }

    return TRUE;
}

//Verify if A is diagonal

int is_diagonal(Matrix *A){

    for(int i = 0; i < A->cols; i++){

        for (int j = 0; j < A->rows; j++) {

            if ( i == j ) {

                if( fabs( A->data[ A->cols*i + j ] ) < 0.000000001 )
                    return FALSE;
            }
        }
    }

    return TRUE;
}

//Verify if A is diagonal

int is_diagonal2(Matrix *A){

    for(int i = 0; i < A->cols; i++){

        for (int j = 0; j < A->rows; j++) {

            if ( i != j ) {

                if( fabs( A->data[ A->cols*i + j ] ) > 0.000000001 )
                    return FALSE;
            }
        }
    }

    return TRUE;
}

Matrix * transpose(Matrix *A){

    // Z = A^T
    Matrix *Z = malloc( sizeof( Z ) );

    // Initialize Matrix struct
    Z->rows = A->cols;
    Z->cols = A->rows;
    Z->data = malloc( Z->rows*Z->cols*sizeof( Z->data ) );

    for(int i = 0; i < Z->rows; i++){

        for(int j = 0; j < Z->cols; j++){

            Z->data[ Z->cols*i + j ] = A->data[ A->cols*j + i ];
        }
    }

    return Z;
}
