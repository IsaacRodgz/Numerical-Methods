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
#include "solve_matrix_direct.h"
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

void powerSolver(Matrix * A, Matrix * eigenVectOld, double* lambdaInit, int num_iters, double epsilon){

    // Flag to indicate Convergence
    int converged = FALSE;

    // Size of matrix and vectors
    int rows = A->rows;

    // Variable to update dominant eigenvalue
    double lambdaNew;

    // Helper variable for dot product of vectors
    double accum;

    // Variable to iterate algorithm
    int i = 0;

    // Vector to calculate v_k
    Matrix *eigenVectNew = malloc( sizeof( eigenVectNew ) );

    // Normalize initial vector v_0
    double norm = vectNorm(eigenVectOld);
    for (int k = 0; k < rows; k++)
        eigenVectOld->data[k] *= (1/norm);

    // Start iterations

    for (i = 0; i < num_iters; i++) {

        // Calculate w =  A * v_(k-1)
        eigenVectNew = multiply(A, eigenVectOld);

        // Calculate v_(k) = w/norm(w)
        norm = vectNorm(eigenVectNew);
        for (int k = 0; k < rows; k++)
            eigenVectOld->data[k] = eigenVectNew->data[k] * (1/norm);

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

        if ( fabs( (*lambdaInit) - lambdaNew ) < epsilon ) {
            (*lambdaInit) = lambdaNew;
            converged = TRUE;
            break;
        }

        (*lambdaInit) = lambdaNew;
    }

    if (converged == TRUE) {
        printf("\nConverged after %d iterations\n\n", i);
    }
    else
        printf("\nMethod did not converge in given iterations. Returning last solution.\n\n");
}

void inversePowerSolver(Matrix * A, Matrix * eigenVect, double* lambdaInit, int num_iters, double epsilon){

    // Flag to indicate Convergence
    int converged = FALSE;

    // Size of matrix and vectors
    int rows = A->rows;

    // Variable to update dominant eigenvalue
    double lambdaNew;

    // Helper variables for dot product of vectors
    double accum;

    // Variable to iterate algorithm
    int i = 0;

    // Factor matrix A
    factor_cholesky_modified(A);

    // Vector to calculate v_(k-1)
    Matrix *eigenVectOld = malloc( sizeof( eigenVectOld ) );
    eigenVectOld->rows = A->rows;
    eigenVectOld->cols = 1;
    eigenVectOld->data = malloc( eigenVectOld->rows * eigenVectOld->cols * sizeof( eigenVectOld->data ) );

    for (int i = 0; i < eigenVectOld->rows; i++) {
        eigenVectOld->data[i] = (double)rand()/RAND_MAX*2.0-1.0;
    }

    // Vector to calculate v_k
    Matrix *eigenVectNew = malloc( sizeof( eigenVectNew ) );

    // Start iterations

    for (i = 0; i < num_iters; i++) {

        // Calculate v_(k) = w/norm(w)
        double norm = vectNorm(eigenVectOld);
        for (int k = 0; k < rows; k++)
            eigenVectOld->data[k] = eigenVectOld->data[k] * (1/norm);

        // Calculate w =  A * v_(k-1)
        eigenVectNew = solve_cholesky_modified(A, eigenVectOld, 0);

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

        if ( fabs( (*lambdaInit) - lambdaNew ) < epsilon ) {
            (*lambdaInit) = lambdaNew;
            converged = TRUE;
            break;
        }

        swap(&eigenVectNew, &eigenVectOld);
        (*lambdaInit) = lambdaNew;
    }

    if (converged == TRUE) {
        printf("----------------------------------------------\n");
        printf("\nConverged after %d iterations\n\n", i);
        for (int i = 0; i < eigenVect->rows; i++) {
            eigenVect->data[i] = eigenVectOld->data[i];
        }
    }
    else
        printf("\nMethod did not converge in given iterations. Returning last solution.\n\n");
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

    factor_cholesky_modified(A);

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
            eigenVectNew = solve_cholesky_modified(A, eigenVectOld, 0);

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

void maxOffDiagonal(Matrix * A, int* p, int* q){

    double max = 0;

    for (int i = 0; i < A->rows; i++) {

        for (int j = 0; j < A->cols; j++) {

            if ( i != j) {

                if( fabs( A->data[ i*A->cols + j] ) > max){

                    max = fabs( A->data[ i*A->cols + j] );
                    p[0] = i;
                    q[0] = j;

                }
            }
        }
    }
}

void jacobiSolver(Matrix * A, Matrix * F, int num_iters, double epsilon){

    for (int i = 0; i < F->rows; i++) {
        for (int j = 0; j < F->cols; j++) {
            if (i == j)
                F->data[i*F->cols + j] = 1.0;
            else
                F->data[i*F->cols + j] = 0.0;
        }
    }

    double* FP = malloc( F->rows * sizeof *FP );
    double* FQ = malloc( F->rows * sizeof *FQ );

    for (int i = 0; i < num_iters; i++) {

        // Get row, column position of max element of A off the diagonal
        int p;
        int q;

        maxOffDiagonal(A, &p, &q);

        if ( fabs( A->data[p*A->cols + q] ) < epsilon ) {
            printf("\nConverged at iteration: %d\n\n", i+1);
            break;
        }

        // Copy columns p and q of matrix of eigenvectors F
        for (int k = 0; k < F->rows; k++) {
            FP[k] = F->data[k*F->cols + p];
            FQ[k] = F->data[k*F->cols + q];
        }

        // Calculate 2*a_ij / ( a_jj - a_ii )
        double w = (A->data[q*(A->cols+1)] - A->data[p*(A->cols+1)]) / (2*A->data[p*A->cols + q]);
        double t = (1/(fabs(w) + sqrt(w*w + 1)))*( w > 0? 1 : -1 );

        double c = 1/sqrt(t*t + 1);
        double s = t*c;
        //double tau = s/(1+c);

        // Update matrix of eigenvectors F
        for (int k = 0; k < F->rows; k++) {
            F->data[k*F->cols + p] = c*FP[k] - s*FQ[k];
            F->data[k*F->cols + q] = s*FP[k] + c*FQ[k];
        }

        double temp_pq = A->data[p*A->cols + q];
        A->data[p*A->cols + q] = 0;
        A->data[q*A->cols + p] = 0;

        double temp_pp = A->data[p*(A->cols+1)];
        A->data[p*(A->cols+1)] = (c*c)*temp_pp + (s*s)*A->data[q*(A->cols+1)] - (2*c*s)*temp_pq;
        A->data[q*(A->cols+1)] = (s*s)*temp_pp + (c*c)*A->data[q*(A->cols+1)] + (2*c*s)*temp_pq;

        for (int j = 0; j < A->rows; j++) {
            double temp;
            if ( j != p & j != q ) {
                temp = A->data[j*A->cols + p];
                A->data[j*A->cols + p] = c*temp - s*A->data[j*A->cols + q];
                A->data[j*A->cols + q] = c*A->data[j*A->cols + q] + s*temp;

                A->data[p*A->cols + j] = A->data[j*A->cols + p];
                A->data[q*A->cols + j] = A->data[j*A->cols + q];

                //temp = A->data[p*A->cols + j];
                //A->data[p*A->cols + j] = c*temp + s*A->data[q*A->cols + j];
                //A->data[q*A->cols + j] = c*A->data[q*A->cols + j] - s*temp;
            }
        }
    }
}
