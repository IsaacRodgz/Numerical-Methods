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

        accum = 0;

        // Calculate A*v_k
        eigenVectNew = multiply(A, eigenVectOld);

        // v_k * (A*v_k)
        for (int j = 0; j < rows; j++) {
            accum += eigenVectNew->data[j] * eigenVectOld->data[j];
        }

        lambdaNew = accum;

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

void inversePowerSolver(Matrix * A, Matrix * eigenVectOld, double* lambdaInit, int num_iters, double epsilon){

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

    int* pivots = malloc( A->rows * sizeof *pivots );
    factor_doolittle_pivoting(A, pivots);

    // Vector to calculate v_k
    Matrix *eigenVectNew = malloc( sizeof( eigenVectNew ) );

    // Normalize initial vector v_0
    double norm = vectNorm(eigenVectOld);
    for (int k = 0; k < rows; k++)
        eigenVectOld->data[k] *= (1/norm);

    // Start iterations

    for (i = 0; i < num_iters; i++) {

        // Calculate w =  A * v_(k-1)
        eigenVectNew = solve_doolittle_pivoting(A, eigenVectOld, pivots);

        // Calculate v_(k) = w/norm(w)
        norm = vectNorm(eigenVectNew);
        for (int k = 0; k < rows; k++)
            eigenVectOld->data[k] = eigenVectNew->data[k] * (1/norm);

        // Calculate A*v_k
        eigenVectNew = multiply(A, eigenVectOld);

        // Calculate dominant eigenvalue and store in lambdaNew

        accum = 0;

        // v_k * (A*v_k)
        for (int j = 0; j < rows; j++) {
            accum += eigenVectNew->data[j] * eigenVectOld->data[j];
        }

        lambdaNew = accum;

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

void kPowerSolver(Matrix * A, Matrix * eigenVectOld, double* lambdaInit, int num_iters, double epsilon){

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

        accum = 0;

        // Calculate A*v_k
        eigenVectNew = multiply(A, eigenVectOld);

        // v_k * (A*v_k)
        for (int j = 0; j < rows; j++) {
            accum += eigenVectNew->data[j] * eigenVectOld->data[j];
        }

        lambdaNew = accum;

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
