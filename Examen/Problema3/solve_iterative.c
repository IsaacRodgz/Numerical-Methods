/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-26
Description: This program includes functions that solve systems A*x = b of n equations with n variables,
through two different iterative methodos methods, namely, Jacobi and Gauss-Seidel.
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "solve_iterative.h"
#include "matrix_struct.h"
#define TRUE 1
#define FALSE 0


double vectNorm(Matrix* vect){

    double sum = 0;

    for (int i = 0; i < vect->rows; i++) {
        sum += ( vect->data[i] * vect->data[i] );
    }

    if ( sum == 0.0 ) {
        fprintf(stderr, "\n[Error] Vector of norm 0 found\n\n");
    }

    return sqrt(sum);

}

void swap(Matrix ** A, Matrix ** B){

    Matrix* temp = *A;
    *A = *B;
    *B = temp;

}

// Helper function to verify if A is symmetric

int is_simetric(Matrix *A){

    for(int i = 0; i < A->rows; i++){

        for(int j = 0; j < A->cols; j++){

            if(A->data[ i*A->rows + j ] != A->data[ j*A->rows + i ]){
                printf("%d, %d\n", i, j);
                printf("%e, %e\n", A->data[ i*A->rows + j ], A->data[ j*A->rows + i ]);
                return FALSE;
            }
        }
    }

    return TRUE;
}

void deflation(Matrix * eigenVects, Matrix * eigenVectInit, int currCol){

    double* temp = calloc( eigenVectInit->rows, sizeof *temp );

    for (int i = 0; i < currCol; i++) {

        double a = 0;

        for (int m = 0; m < eigenVects->cols; m++) {

            a += eigenVectInit->data[m] * eigenVects->data[ i*eigenVects->cols + m ];
        }

        for (int m = 0; m < eigenVects->cols; m++) {

            temp[m] += a * eigenVects->data[ i*eigenVects->cols + m ];
        }
    }

    for (int i = 0; i < eigenVectInit->rows; i++) {

        eigenVectInit->data[i] -= temp[i];
    }
}

void kPowerSolver(Matrix * A, Matrix * eigenVects, Matrix * eigenVals, int num_iters, double epsilon, int k){

    if ( is_simetric(A) == FALSE ) {
        printf("Matrix is not symmetric, cannot apply algorithm\n");
        exit(-1);
    }

    // Size of matrix and vectors
    int rows = A->rows;

    // Variable to update dominant eigenvalue
    double lambdaNew;

    // Vector to calculate v_k
    Matrix *eigenVectNew = malloc( sizeof( eigenVectNew ) );

    // Vector to calculate v_(k-1)
    Matrix *eigenVectOld = malloc( sizeof( eigenVectOld ) );
    eigenVectOld->rows = A->rows;
    eigenVectOld->cols = 1;
    eigenVectOld->data = malloc( eigenVectOld->rows * eigenVectOld->cols * sizeof( eigenVectOld->data ) );

    for (int s = 0; s < k; s++) {

        // Flag to indicate Convergence
        int converged = FALSE;

        // Initialize and Normalize initial vector v_0

        for (int i = 0; i < eigenVectOld->rows; i++) {
            eigenVectOld->data[i] = (double)rand()/RAND_MAX*2.0-1.0;
        }

        double norm = vectNorm(eigenVectOld);
        for (int k = 0; k < rows; k++)
            eigenVectOld->data[k] *= (1/norm);

        // Initialize lambda_0
        eigenVals->data[s] = 0;

        // Start iterations
        int i = 0;
        for (i = 0; i < num_iters; i++) {

            if (s>0)
                deflation(eigenVects, eigenVectOld, s);

            // Calculate w =  A * v_(k-1)
            eigenVectNew = multiply(A, eigenVectOld);

            // Calculate v_(k) = w/norm(w)
            norm = vectNorm(eigenVectNew);
            for (int k = 0; k < rows; k++)
                eigenVectOld->data[k] = eigenVectNew->data[k] * (1/norm);

            // Calculate dominant eigenvalue and store in lambdaNew

            // Calculate A*v_k
            eigenVectNew = multiply(A, eigenVectOld);

            double accum = 0;

            // v_k * (A*v_k)
            for (int j = 0; j < rows; j++) {
                accum += eigenVectNew->data[j] * eigenVectOld->data[j];
            }

            lambdaNew = accum;

            // Check for convergence and stop or update the eigenvalue

            if ( fabs( eigenVals->data[s] - lambdaNew ) < epsilon ) {
                eigenVals->data[s] = lambdaNew;
                converged = TRUE;
                break;
            }

            eigenVals->data[s] = lambdaNew;
        }

        if (converged == TRUE) {
            printf("----------------------------------------------\n");
            printf("\nConverged for eigenvector_%d after %d iterations\n\n", s+1, i+1);

            for (int i = 0; i < eigenVects->cols; i++) {
                eigenVects->data[ s*eigenVects->cols + i ] = eigenVectOld->data[i];
            }
        }
        else
            printf("\nMethod did not converge in given iterations. Returning last solution.\n\n");
    }
}

void kInversePowerSolver(Matrix * A, Matrix * eigenVects, Matrix * eigenVals, int num_iters, double epsilon, int k){

    if ( is_simetric(A) == FALSE ) {
        printf("Matrix is not symmetric, cannot apply algorithm\n");
        exit(-1);
    }

    // Size of matrix and vectors
    int rows = A->rows;

    // Variable to update dominant eigenvalue
    double lambdaNew;

    // Vector to calculate v_k
    Matrix *eigenVectNew = malloc( sizeof( eigenVectNew ) );

    // Vector to calculate v_(k-1)
    Matrix *eigenVectOld = malloc( sizeof( eigenVectOld ) );
    eigenVectOld->rows = A->rows;
    eigenVectOld->cols = 1;
    eigenVectOld->data = malloc( eigenVectOld->rows * eigenVectOld->cols * sizeof( eigenVectOld->data ) );

    int* pivots = malloc(A->rows * sizeof *pivots);
    factor_doolittle_pivoting(A, pivots);

    for (int s = 0; s < k; s++) {

        // Flag to indicate Convergence
        int converged = FALSE;

        // Initialize and Normalize initial vector v_0

        for (int i = 0; i < eigenVectOld->rows; i++) {
            eigenVectOld->data[i] = (double)rand()/RAND_MAX*2.0-1.0;
        }

        double norm;

        // Initialize lambda_0
        eigenVals->data[s] = 0;

        // Start iterations
        int i = 0;
        for (i = 0; i < num_iters; i++) {

            if (s>0)
                deflation(eigenVects, eigenVectOld, s);

            // Calculate v_(k) = w/norm(w)
            norm = vectNorm(eigenVectOld);
            for (int k = 0; k < rows; k++)
                eigenVectOld->data[k] = eigenVectOld->data[k] * (1/norm);

            // Calculate w =  A * v_(k-1)
            eigenVectNew = solve_doolittle_pivoting(A, eigenVectOld, pivots);

            // Calculate dominant eigenvalue and store in lambdaNew

            double accum2 = 0;
            double accum = 0;

            // v_k * (A*v_k)
            for (int j = 0; j < rows; j++) {
                accum += eigenVectNew->data[j] * eigenVectOld->data[j];
                accum2 += eigenVectNew->data[j] * eigenVectNew->data[j];
            }

            lambdaNew = accum/accum2;

            // Check for convergence and stop or update the eigenvalue

            if ( fabs( eigenVals->data[s] - lambdaNew ) < epsilon ) {
                eigenVals->data[s] = lambdaNew;
                converged = TRUE;
                break;
            }

            eigenVals->data[s] = lambdaNew;
            swap(&eigenVectNew, &eigenVectOld);
        }

        if (converged == TRUE) {
            printf("----------------------------------------------\n");
            printf("\nConverged for eigenvector_%d after %d iterations\n\n", s+1, i+1);

            for (int i = 0; i < eigenVects->cols; i++) {
                eigenVects->data[ s*eigenVects->cols + i ] = eigenVectOld->data[i];
            }
        }
        else
            printf("\nMethod did not converge in given iterations. Returning last solution.\n\n");
    }
}

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

