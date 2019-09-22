/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-09-18
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "solve_matrix_direct.h"
#include "solve_iterative3.h"
#include "solve_iterative2.h"
#include "matrix_struct.h"
#define TRUE 1
#define FALSE 0

void rayleighSolver(Matrix * A, Matrix * eigenVect, double * lambdaInit, int num_iters, double epsilon){

    // Make a copy of A
    Matrix *A_shifted = malloc( sizeof( A_shifted ) );
    A_shifted->rows = A->rows;
    A_shifted->cols = A->cols;
    A_shifted->data = malloc( A_shifted->rows * A_shifted->cols * sizeof( A_shifted->data ) );

    copy(A, A_shifted);

    // Calculate A - sigma*I. Shifted A.
    for (int i = 0; i < A_shifted->rows; i++) {
        A_shifted->data[(A_shifted->rows+1)*i] -= *(lambdaInit);
    }

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
    int* pivots = malloc(A_shifted->rows * sizeof *pivots);
    factor_doolittle_pivoting(A_shifted, pivots);

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

        // Calculate w =  A * v_(k-1)
        eigenVectNew = solve_doolittle_pivoting(A_shifted, eigenVectOld, pivots);

        // Calculate v_(k) = w/norm(w)
        double norm = vectNorm(eigenVectNew);
        for (int k = 0; k < rows; k++)
            eigenVectOld->data[k] = eigenVectNew->data[k] * (1/norm);

        // Calculate dominant eigenvalue and store in lambdaNew

        #pragma omp parallel for
        for (int j = 0; j < A->rows; ++j){

            double sum = 0;

            for (int k = 0; k < A->rows; ++k){

                sum += A->data[ j*A->cols + k ]*eigenVectOld->data[k];
            }
            eigenVectNew->data[j] = sum;
        }

        accum = 0;

        #pragma omp parallel for reduction(+:accum)
        for (int j = 0; j < eigenVectNew->rows; ++j)
            accum += eigenVectNew->data[j] * eigenVectOld->data[j];

        lambdaNew = accum;

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
        //printf("----------------------------------------------\n");
        printf("\nConverged after %d iterations\n\n", i);
        for (int i = 0; i < eigenVect->rows; i++) {
            eigenVect->data[i] = eigenVectOld->data[i];
        }
    }
    else
        printf("\nMethod did not converge in given iterations. Returning last solution.\n\n");

    free(A_shifted->data);
    free(A_shifted);
    free(eigenVectOld->data);
    free(eigenVectOld);
    free(eigenVectNew->data);
    free(eigenVectNew);
}

void subspaceSolver(Matrix * A, Matrix * FI, Matrix * LA, int num_iters, double epsilon, int k) {

    // FI[n][p] - matrix of eigenvectors of A
    // LA[p][p] - matrix of eigenvalues of A

    // Q[n][p]

    Matrix *Q = malloc( sizeof( Q ) );
    Q->rows = FI->rows;
    Q->cols = FI->cols;
    Q->data = malloc( Q->rows * Q->cols * sizeof( Q->data ) );

    // R[p][p]

    Matrix *R = malloc( sizeof( R ) );
    R->rows = LA->rows;
    R->cols = LA->cols;
    R->data = malloc( R->rows * R->cols * sizeof( R->data ) );

    // Temp matrix for multiply( Q^T, A )

    Matrix *T = malloc( sizeof( T ) );
    T->rows = Q->cols;
    T->cols = A->cols;
    T->data = malloc( T->rows * T->cols * sizeof( T->data ) );

    // Initialize matrix of eigenvectors FI

    for (size_t i = 0; i < FI->rows; i++) {
        for (int j = 0; j < FI->cols; j++) {
            FI->data[FI->cols*i + j] = (double)rand()/RAND_MAX*2.0-1.0;
        }
    }

    for (int i = 0; i < num_iters; i++) {

        // Iterate in subspace

        QRFactor(FI, Q, R);

        // T = multiply( Q^T, A )

        int l, k;
        #pragma omp parallel for private(l,k)
        for(int j = 0; j < T->rows; j++){
            for(k = 0; k < T->cols; k++){

                T->data[ T->cols*j + k ] = 0;

                for(l = 0; l < A->rows; l++){
                    T->data[ T->cols*j + k ] += Q->data[ Q->cols*l + j ] * A->data[ A->cols*l + k ];
                }
            }
        }

        // LA = multiply( T, Q );

        #pragma omp parallel for private(l,k)
        for(int j = 0; j < LA->rows; j++){
            for(k = 0; k < LA->cols; k++){

                LA->data[ LA->cols*j + k ] = 0;

                for(l = 0; l < Q->rows; l++){
                    LA->data[ LA->cols*j + k ] += T->data[ T->cols*j + l ] * Q->data[ Q->cols*l + k ];
                }
            }
        }

        // FI = multiply( A, Q );

        #pragma omp parallel for private(l,k)
        for(int j = 0; j < FI->rows; j++){
            for(k = 0; k < FI->cols; k++){

                FI->data[ FI->cols*j + k ] = 0;

                for(l = 0; l < Q->rows; l++){
                    FI->data[ FI->cols*j + k ] += A->data[ A->cols*j + l ] * Q->data[ Q->cols*l + k ];
                }
            }
        }


        if (is_diagonal2(LA) == TRUE) {
            printf("\nSubspace method converged after %d iterations\n", i);
            return;
        }
    }

    printf("\nSubspace method did not converged after %d iterations. Returning last computations\n", num_iters);
}

void cGradientSolver(Matrix * A, Matrix * b, Matrix * x, int numIters, double epsilon){

    // Verify that A is simetric

    if( is_simetric(A) == FALSE){

        printf("\nA is not symmetric. System cannot be solved by conjugate gradient.\n\n");
        exit(-1);
    }

    // Temp array for p_k.T * A * p_k

    double* temp = malloc( A->rows * sizeof *temp );

    // Calculate A * x

    #pragma omp parallel for
    for (int j = 0; j < A->rows; ++j){

        double sum = 0;

        for (int k = 0; k < A->rows; ++k){

            sum += A->data[ j*A->cols + k ]*x->data[k];
        }
        temp[j] = sum;
    }

    Matrix *r = malloc( sizeof( r ) );
    r->rows = b->rows;
    r->cols = b->cols;
    r->data = malloc( r->rows * r->cols * sizeof( r->data ) );

    // Initialize r_0 = b - A*x

    for (int i = 0; i < r->rows; i++) {
        r->data[i] = b->data[i] - temp[i];
    }

    Matrix *p = malloc( sizeof( p ) );
    p->rows = r->rows;
    p->cols = r->cols;
    p->data = malloc( p->rows * p->cols * sizeof( p->data ) );

    // Initialize p_0 = r_0

    for (int i = 0; i < p->rows; i++) {
        p->data[i] = r->data[i];
    }

    double alpha = 0;
    double beta = 0;

    double rNorm_old = 0;

    #pragma omp parallel for reduction(+:rNorm_old)
    for (int j = 0; j < r->rows; ++j)
        rNorm_old += r->data[j] * r->data[j];

    for (int i = 0; i < numIters; i++) {

        // Calculate A * p_k

        #pragma omp parallel for
        for (int j = 0; j < A->rows; ++j){

            double sum = 0;

            for (int k = 0; k < A->rows; ++k){

                sum += A->data[ j*A->cols + k ]*p->data[k];
            }
            temp[j] = sum;
        }

        // Calculate alpha_k = (p_k * r_k) / (p_k * A * p_k)

        double denom = 0;
        #pragma omp parallel for reduction(+:denom)
        for (int j = 0; j < p->rows; ++j)
            denom += p->data[j] * temp[j];

        alpha = rNorm_old/denom;

        // Update x_k+1 = x_k + alpha*p_k

        for (int j = 0; j < x->rows; j++) {
            x->data[j] += alpha*p->data[j];
        }

        double rNorm_new = 0;

        // Update r_k+1 = r_k - alpha*A*p_k
        for (int j = 0; j < r->rows; j++) {
            r->data[j] -= alpha*temp[j];
            rNorm_new += r->data[j] * r->data[j];
        }

        if ( sqrt(rNorm_new) < epsilon ) {
            printf("\nAlgorithm converged after %d iterations\n\n", i+1);
            return;
        }

        beta = rNorm_new/rNorm_old;

        // Update p_k+1 = r_k+1 + beta*p_k

        for (int j = 0; j < p->rows; j++) {
            p->data[j] = r->data[j] + beta*p->data[j];
        }

        rNorm_old = rNorm_new;
    }

    printf("\nAlgorithm could not converge after %d iterations\n", numIters);
}

void QRFactor(Matrix * A, Matrix * Q, Matrix * R){

    // Calculate r_00

    double norm = 0;

    for (int i = 0; i < A->rows; i++) {
        norm += A->data[A->cols*i] * A->data[A->cols*i];
    }

    norm = sqrt(norm);
    R->data[0] = norm;

    // Calculate q_0

    for (int i = 0; i < Q->rows; i++) {
        Q->data[Q->cols*i] = A->data[A->cols*i] * (1/norm);
    }

    // Vector to form new q's
    double* a_temp = malloc( A->rows * sizeof *a_temp );

    for (int i = 1; i < A->cols; i++) {

        // Initialize a_temp to column A[i]
        for (int j = 0; j < A->rows; j++) {
            a_temp[j] = A->data[ A->cols*j + i ];
        }

        for (int j = 0; j < i; j++) {

            double dot = 0;

            for (int k = 0; k < A->rows; k++) {
                //printf("Q[%d], A[%d]\n", Q->cols*k + j, A->cols*k + size);
                dot += Q->data[ Q->cols*k + j ] * A->data[ A->cols*k + i ];
            }

            R->data[ Q->cols*j + i ] = dot;

            for (int k = 0; k < A->rows; k++) {
                a_temp[k] -= dot * Q->data[ Q->cols*k + j];
            }
        }

        norm = 0;
        for (int j = 0; j < A->rows; j++) {
            norm += a_temp[j] * a_temp[j];
        }
        norm = sqrt(norm);

        // Calculate q_i

        for (int j = 0; j < Q->rows; j++) {
            Q->data[Q->cols*j + i] = a_temp[j] * (1/norm);
        }

        // Calculate r_ii

        R->data[ (Q->cols+1)*i ] = 0;

        for (int j = 0; j < A->rows; j++) {
            //printf("Q[%d], A[%d]\n", Q->cols*k + j, A->cols*k + size);
            R->data[ (Q->cols+1)*i ] += Q->data[ Q->cols*j + i ] * A->data[ A->cols*j + i ];
        }
    }

    free(a_temp);
}

void QRSolve(Matrix * A,  Matrix * FI, int numIters, double epsilon){

    Matrix *Q = malloc( sizeof( Q ) );
    Q->rows = A->rows;
    Q->cols = A->cols;
    Q->data = malloc( Q->rows * Q->cols * sizeof( Q->data ) );

    Matrix *R = malloc( sizeof( R ) );
    R->rows = A->rows;
    R->cols = A->cols;
    R->data = malloc( R->rows * R->cols * sizeof( R->data ) );

    // Initialize FI with Identity matrix

    for (int i = 0; i < FI->rows; i++) {
        for (int j = 0; j < FI->cols; j++) {
            if(i!=j)
                FI->data[i*FI->cols + j] = 0.0;
            else
                FI->data[i*FI->cols + j] = 1.0;
        }
    }

    for (int i = 0; i < numIters; i++) {

        QRFactor(A, Q, R);

        // Update eigenvalues

        int j, k, l;
        #pragma omp parallel for private(k,l)
        for(j = 0; j < A->rows; j++){
            for(k = 0; k < A->cols; k++){

                A->data[ A->cols*j + k ] = 0;

                for(l = 0; l < Q->rows; l++){
                    A->data[ A->cols*j + k ] += R->data[ R->cols*j + l ] * Q->data[ Q->cols*l + k ];
                }
            }
        }

        // Update eigenvectors

        #pragma omp parallel for private(k,l)
        for(j = 0; j < FI->rows; j++){
            for(k = 0; k < FI->cols; k++){

                FI->data[ FI->cols*j + k ] = 0;

                for(l = 0; l < Q->rows; l++){
                    FI->data[ FI->cols*j + k ] += FI->data[ FI->cols*j + l ] * Q->data[ Q->cols*l + k ];
                }
            }
        }
    }
}
