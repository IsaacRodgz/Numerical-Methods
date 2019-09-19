/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-09-18
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
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

void cGradientSolver(Matrix * A, Matrix * b, Matrix * x, int numIters, double epsilon){

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

        // Calculate beta_k = () / ()
        /*
        num = 0;
        denom = 0;

        #pragma omp parallel for reduction(+:num)
        for (int j = 0; j < r->rows; ++j)
            num += r->data[j] * r->data[j];

        #pragma omp parallel for reduction(+:denom)
        for (int j = 0; j < p->rows; ++j)
            denom += p->data[j] * p->data[j];
        */
        beta = rNorm_new/rNorm_old;

        // Update p_k+1 = r_k+1 + beta*p_k

        for (int j = 0; j < p->rows; j++) {
            p->data[j] = r->data[j] + beta*p->data[j];
        }

        rNorm_old = rNorm_new;
    }

    printf("\nAlgorithm could not converge after %d iterations\n", numIters);
}

void QRFactor(Matrix * A, Matrix * Q, Matrix * R, int numIters, double epsilon){
    /*
    for (int i = 0; i < A->rows; i++) {

        double* v = malloc( (A->rows-i) * sizeof *v );

        int curr_row = i;

        for (int j = 0; j < (A->rows-i); j++) {
            v[j] = A->data[ A->rows*curr_row + i ];
            curr_row++;
        }

        double norm = 0;
        for (int j = 0; j < (A->rows-i); j++) {
            norm += v[j]*v[j];
        }
        norm = sqrt(norm);

        v[0] += norm*( v[0]>=0? 1 : -1 );

        free(v);
    }
    */

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

    // Find remaining of R and Q
    int size = 1;

    // Vector to form new q's
    double* a_temp = malloc( A->rows * sizeof *a_temp );

    for (int i = 1; i < R->cols; i++) {

        // Initialize a_temp to column A[i]
        for (int j = 0; j < A->rows; j++) {
            a_temp[j] = A->data[ A->cols*j + size ];
        }

        // Sum{q*r}
        double sum = 0;

        for (int j = 0; j < size; j++) {

            double dot = 0;

            for (int k = 0; k < A->rows; k++) {
                //printf("Q[%d], A[%d]\n", Q->cols*k + j, A->cols*k + size);
                dot += Q->data[ Q->cols*k + j ] * A->data[ A->cols*k + size ];
            }

            R->data[ Q->cols*j + size ] = dot;

            for (int k = 0; k < A->rows; k++) {
                a_temp[j] -= dot * Q->data[ Q->cols*k + j];
            }

            //printf("--\n");
        }
        //printf("\n\n");

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

        R->data[ (Q->cols+1)*size ] = norm;

        size++;
    }
}
