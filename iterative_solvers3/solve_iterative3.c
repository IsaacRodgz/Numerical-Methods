/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-26
Description: This program includes functions that solve systems A*x = b of n equations with n variables,
through two different iterative methodos methods, namely, Jacobi and Gauss-Seidel.
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "solve_iterative3.h"
#include "solve_iterative2.h"
#include "matrix_struct.h"
#define TRUE 1
#define FALSE 0

void rayleighSolver(Matrix * A, Matrix * eigenVect, double * lambda, int num_iters, double epsilon){

    Matrix *r;

    Matrix *Ap = malloc( sizeof( Ap ) );
    Ap->rows = A->rows;
    Ap->cols = A->cols;
    Ap->data = malloc( Ap->rows * Ap->cols * sizeof( Ap->data ) );

    copy(A, Ap);

    double norm = vectNorm(eigenVect);
    for (int k = 0; k < eigenVect->rows; k++) {
        eigenVect->data[k] *= eigenVect->data[k]*(1/norm);
    }

    int i;
    for (i = 0; i < num_iters-1; i++) {

        for (int k = 0; k < Ap->cols; k++) {
            Ap->data[k*(Ap->cols+1)] = A->data[k*(A->cols+1)] - (*lambda);
        }

        r = multiply(Ap, eigenVect);

        double dot_prod = 0;
        for (int k = 0; k < eigenVect->rows; k++) {
            dot_prod += eigenVect->data[k] * r->data[k];
        }

        if ( fabs(dot_prod) < epsilon ) {
            break;
        }

        (*lambda) += dot_prod;

    }

    printf("\nSe llego a la convergencia en %d iteraciones\n", i+1);

    free(r->data);
    free(r);
    free(Ap->data);
    free(Ap);
}

void subspaceSolver(Matrix * A, Matrix * FI, Matrix * LA, int num_iters, double epsilon, int k) {

    // FI - matrix of eigenvectors of A
    // LA - matrix of eigenvalues of A

    // Get initial guess with power method and deflation
    printf("\n\n**************************\n");
    printf("Solving by power iteration");
    printf("\n\n**************************\n");
    kPowerSolver(A, FI, LA, num_iters, epsilon, k);

    // B[m][m] = FI[m][n]*A[n][n]*FI[n][m]

    Matrix *B = malloc( sizeof( B ) );
    B->rows = FI->rows;
    B->cols = FI->rows;
    B->data = malloc( B->rows * B->cols * sizeof( B->data ) );

    // Factor B by jacobi: B[m][m] -> Q[m][m]*L[m][m]*Q[m][m]

    Matrix *Q = malloc( sizeof( Q ) );
    Q->rows = B->rows;
    Q->cols = B->cols;
    Q->data = malloc( Q->rows * Q->cols * sizeof( Q->data ) );

    // Iterate with Jacobi to update Eigenvectors FI with FI[m][n] = Q[m][m]*FI[m][n]
    printf("\n*********************\n");
    printf("Iterating in subspace");
    printf("\n*********************\n\n");
    for (int i = 0; i < num_iters; i++) {

        B = multiply( multiply( FI, A ), transpose(FI) );

        if (is_diagonal(B) == TRUE) {
            printf("\nMethod converged after %d iterations\n", i);
            break;
        }

        printf("\n*********************\n");
        printf("Iterating Jacobi");
        printf("\n*********************\n\n");
        jacobiSolver(B, Q, num_iters, epsilon);

        FI = multiply(Q, FI);
    }

}
