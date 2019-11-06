/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-12
Description: Pogram that solves system A*x = b of n equations with n variables,
through different methods.
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "solve_matrix_direct.h"
#include "matrix_struct.h"
#define TRUE 1
#define FALSE 0

//Verify if A is diagonal

int is_diagonal(Matrix *A){

    int cols = A->rows;

    for(int i = 0; i < cols; i++){

        if( fabs( A->data[ (cols + 1)*i ] ) == 0.0 )
            return FALSE;
    }

    return TRUE;
}

// Solve linear system where A is a lower triangular matrix

Matrix * solve_lower_triang(Matrix *A, Matrix *b, int fill_diag){

    if( A->rows != A->cols ){

        printf("\nNumber of rows and columns don't match. System cannot be solved\n\n");
        exit(-1);
    }

    if( is_diagonal(A) == FALSE ){

        printf("\nThere are diagonal elements equal to zero. System cannot be solved.\n\n");
        exit(-1);
    }

    double sum;
    int rows = A->rows;
    int cols = A->cols;

    // Vector of solution
    Matrix *x = malloc( sizeof( x ) );

    // Initialize Matrix struct
    x->rows = rows;
    x->cols = 1;
    x->data = malloc( x->rows*x->cols*sizeof( x->data ) );

    // fill_diag = 1 changes diagional of A (A[i][i]) to 1's

    if(fill_diag == 1)
        x->data[0] = b->data[0];
    else
        x->data[0] = b->data[0]/A->data[0];

    // Calculate x[1] to x[n-1] with formula
    // x[i] = ( b[i] - Sum_(0 to i-1)_(A[i][j]*x[j]) ) / (A[i][i])

    for(int i = 1; i < rows; i++){

        sum = 0;

        for(int j = 0; j < i; j++){

            sum += A->data[ cols*i + j ] * x->data[ j ];
        }

        if(fill_diag == 1)
            x->data[i] = b->data[i] - sum;
        else
            x->data[i] = (b->data[i] - sum) / ( A->data[ i*(rows+1) ] ); // A[i][i] = A->data[ i*(rows+1) ]
    }

    return x;
}

// Solve linear system where A is a upper triangular matrix

Matrix * solve_upper_triang(Matrix *A, Matrix *b, int fill_diag){

    if( A->rows != A->cols ){

        printf("\nNumber of rows and columns don't match. System cannot be solved\n\n");
        exit(-1);
    }

    if( is_diagonal(A) == FALSE ){

        printf("\nThere are diagonal elements equal to zero. System cannot be solved.\n\n");
        exit(-1);
    }

    double sum;
    int rows = A->rows;

    // Vector of solution
    Matrix *x = malloc( sizeof( x ) );

    // Initialize Matrix struct
    x->rows = rows;
    x->cols = 1;
    x->data = malloc( x->rows*x->cols*sizeof( x->data ) );

    // fill_diag = 1 changes diagional of A (A[i][i]) to 1's

    if(fill_diag == 1)
        x->data[rows-1] = b->data[rows-1];
    else
        x->data[rows-1] = b->data[rows-1] / A->data[ rows*rows-1 ];

    // Calculate x[n-2] to x[0] with formula
    // x[i] = ( b[i] - Sum_(i+1 to n-1)_(A[i][j]*x[j]) ) / (A[i][i])

    for(int i = rows-2; i >= 0; i--){

        sum = 0;

        for(int j = i+1; j < rows; j++){

            sum += A->data[ rows*i + j ] * x->data[j];
        }

        if(fill_diag == 1)
            x->data[i] = b->data[i] - sum;
        else
            x->data[i] = (b->data[i] - sum) / A->data[ i*(rows+1) ]; // A[i][i] = A->data[ i*(n+1) ]
    }

    return x;
}

// Helper function to exchange rows from index: i_switch to index: index_set

void set_pivot_row(Matrix *A, int* pivot, int i_switch, int index_set){

    int cols = A->cols;

    // Helper array to switch row and columns of A

    double buffer;

    // Switch rows in A

    for(int i = 0; i < cols; i++){

        buffer = A->data[ i_switch*cols + i ];

        A->data[ i_switch*cols + i ] = A->data[ index_set*cols + i ];

        A->data[ index_set*cols + i ] = buffer;
    }

    int temp = pivot[i_switch];
    pivot[i_switch] = pivot[index_set];
    pivot[index_set] = temp;

}

// Helper function to get the indexes of maximum element in current colum (in absolute value) of matrix A

int max_pivot_index_column(Matrix *A, int limit){

    int cols = A->cols;
    int row = limit;
    double max = 0;

    for(int i = limit; i < cols; i++){

        if( fabs( A->data[ i*cols + limit] ) > max){

            max = fabs( A->data[ i*cols + limit] );
            row = i;

        }
    }

    return row;
}

// Helper function to factor matrix A into L*U with Doolittle and partial pivoting

void factor_doolittle_pivoting(Matrix *A, int *pivot){

    int cols = A->cols;
    int current_pivot;

    for (int i = 0; i < cols; i++)
        pivot[i] = i;

    for (int i = 0; i < cols; i++) {

        current_pivot = max_pivot_index_column(A, i);

        if ( A->data[ current_pivot*(cols+1) ] == 0.0 ) {
            printf("\nThere are diagonal elements equal to zero. System cannot be solved\n\n");
            exit(-1);
        }

        if ( current_pivot != i ) {

            set_pivot_row(A, pivot, current_pivot, i);
        }

        for (int j = i + 1; j < cols; j++) {

            A->data[ j*cols + i ] = A->data[ j*cols + i ] / A->data[ i*(cols+1) ];

            for(int k = i + 1; k < cols; k++){

                A->data[ j*cols + k ] -=  ( A->data[ j*cols + i ] * A->data[ i*cols + k ] );
            }
        }

    }
}

// Helper function to factor matrix A into L*U with Doolittle

void factor_doolittle(Matrix *A){

    int cols = A->cols;

    // Alternate calculation of rows of U and columns of L

    for(int i = 0; i < cols; i++){

        // Calculate row i of U and save in A

        for(int j = i; j < cols; j++){

            // U[i][j] = A[i][j] - Sum_(k=0 to i-1){ L[i][k] * U[k][j] }

            for(int k = 0; k < i; k++){

                // A[i][j] -= A[i][k] * A[k][j]

                A->data[ i*cols + j] -= A->data[ i*cols + k ] * A->data[ k*cols + j ];
            }
        }

        // Calculate column i of L and save in A

        for(int j = i+1; j < cols; j++){

            // L[j][i] = A[j][i]/U[i][i] - Sum_(k=0 to i-1){ L[j][k] * U[k][i] }/U[i][i]

            // A[j][i] = A[j][i]/A[i][i]

            A->data[ j*cols + i ] = A->data[ j*cols + i ] / A->data[ i*(cols+1) ];

            for(int k = 0; k < i; k++){

                // A[j][i] -= (A[j][k] * A[k][i])/A[i][i]

                A->data[ j*cols + i ] -= ( A->data[ j*cols + k] * A->data[ k*cols + i ] ) / A->data[ i*(cols+1) ];
            }
        }
    }
}

Matrix * solve_doolittle_pivoting(Matrix *A, Matrix *b, int *pivots){

    int rows = A->rows;

    // Order b in x_solve

    Matrix *x_solve = malloc( sizeof( x_solve ) );

    x_solve->rows = rows;
    x_solve->cols = 1;
    x_solve->data = malloc( x_solve->rows*x_solve->cols*sizeof( x_solve->data ) );

    for(int i = 0; i < rows; i++)
        x_solve->data[ i ] = b->data[ pivots[i] ];

    // Solve L*y = b where b is x_solve

    x_solve = solve_lower_triang(A, x_solve, 1);

    // Solve U*x = y where y is x_solve

    x_solve = solve_upper_triang(A, x_solve, 0);

    return x_solve;
}

Matrix * solve_doolittle(Matrix *A, Matrix *b, int factor_flag){

    int rows = A->rows;

    // Factor A as LU through Doolittle's method if factor_flag = 1

    if(factor_flag == 1)
        factor_doolittle(A);

    // Solve L*y = b

    Matrix *y = malloc( sizeof( y ) );

    y->rows = rows;
    y->cols = 1;
    y->data = malloc( y->rows*y->cols*sizeof( y->data ) );

    y = solve_lower_triang(A, b, 1);

    // Solve U*x = y

    Matrix *x_solve = malloc( sizeof( x_solve ) );

    x_solve->rows = rows;
    x_solve->cols = 1;
    x_solve->data = malloc( x_solve->rows*x_solve->cols*sizeof( x_solve->data ) );

    x_solve = solve_upper_triang(A, y, 0);

    free(y->data);
    free(y);

    return x_solve;
}

// Find inverse of matrix A through Doolittle

Matrix * solve_inverse(Matrix *A){

    int cols = A->cols;
    int rows = A->rows;

    // Matrix to store inverse of A

    Matrix *inverse_A = malloc( sizeof( inverse_A ) );

    inverse_A->rows = rows;
    inverse_A->cols = cols;
    inverse_A->data = malloc( inverse_A->rows*inverse_A->cols*sizeof( inverse_A->data ) );

    // Array of b. Takes columns of identity matrix

    Matrix *b = malloc( sizeof( b ) );

    b->rows = rows;
    b->cols = 1;
    b->data = malloc( b->rows*b->cols*sizeof( b->data ) );

    // Array to store solution of A*x = b, where x represents each column of A inverse

    Matrix *x = malloc( sizeof( x ) );

    x->rows = rows;
    x->cols = 1;
    x->data = malloc( x->rows*x->cols*sizeof( x->data ) );

    // Factor A to LU through doolittle's algorithm

    factor_doolittle(A);

    // Solve n linear sistems

    for(int i = 0; i < cols; i++){

        // Initialize b to i'th column of identity matrix

        for(int j = 0; j < cols; j++){

            if(j == i)
                b->data[j] = 1;
            else
                b->data[j] = 0;
        }

        // Solve for x with A already factored as LU

        x = solve_doolittle(A, b, 0);

        // Fill i'th column of A inverse with x

        for(int j = 0; j < cols; j++){

            inverse_A->data[ j*cols + i] = x->data[j];
        }
    }

    return inverse_A;
}

