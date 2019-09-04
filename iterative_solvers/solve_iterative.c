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

Matrix * jacobiSolver(Matrix * A, Matrix * b, int num_iters, double tolerance){

    int rows = A->rows;

    // initialize x_new

    Matrix *x_new = malloc( sizeof( x_new ) );

    x_new->rows = rows;
    x_new->cols = 1;
    x_new->data = malloc( x_new->rows * x_new->cols * sizeof( x_new->data ) );

    // initialize x_old

    Matrix *x_old = malloc( sizeof( x_old ) );

    x_old->rows = rows;
    x_old->cols = 1;
    x_old->data = malloc( x_old->rows * x_old->cols * sizeof( x_old->data ) );

    // Helper variable for Sum(j=1 to n){ A[i][j] * x_old[j] }

    double sum;

    // Var to measure convergence between x_new and x_old

    double error;

    // Helper vars to calculate l2-norms form error

    double num_norm;
    double den_norm;

    // Initialize x_old with b

    copy( b, x_old );

    for(int iter = 0; iter < num_iters; iter++){

        error = 0;
        num_norm = 0;
        den_norm = 0;

        for(int i = 0; i < rows; i++){

            // Validate A[i][i] is not zero

            if( A->data[ i*(rows+1) ] == 0 ){

                printf("\nDivision by 0. Cannot continue. Stopped at iteration: %d\n\n", iter+1);
            }

            sum = 0;

            for(int j = 0; j < rows; j++){

                if(j != i)
                    sum += A->data[ i*rows + j ] * x_old->data[j];
            }

            x_new->data[i] = ( b->data[i] - sum ) / A->data[ i*(rows+1) ];

            num_norm += (x_new->data[i] - x_old->data[i]) * (x_new->data[i] - x_old->data[i]);
            den_norm += (x_new->data[i])*(x_new->data[i]);
        }

        error = sqrt(num_norm/den_norm);

        // If convergence is reached, stop iterating and return solution x_new

        if( error < tolerance ){
            printf("\n***********\n\nConvergence reached at iteration %d\n\n***********\n", iter+1);
            return x_new;
        }

        // update x_old with x_new content and continue iteration

        copy(x_new, x_old);
    }

    printf("\nSolution x did not converge\n");

    return x_new;
}


Matrix * gaussSeidelSolver(Matrix * A, Matrix * b, int num_iters, double tolerance){

    int rows = A->rows;

    // initialize x_new

    Matrix *x_new = malloc( sizeof( x_new ) );

    x_new->rows = rows;
    x_new->cols = 1;
    x_new->data = malloc( x_new->rows * x_new->cols * sizeof( x_new->data ) );

    // Helper variable for Sum(j=1 to n){ A[i][j] * x_old[j] }

    double sum;

    // Var to measure convergence between x_new and x_old

    double error;

    // Helper vars to calculate l2-norms form error

    double num_norm;
    double den_norm;

    // Helper variable to calculate error ( x_new[i] - x_old[i] ) without needing X_old

    double temp_x;

    // Initialize x_old with b

    copy( b, x_new );

    for(int iter = 0; iter < num_iters; iter++){

        error = 0;
        num_norm = 0;
        den_norm = 0;

        for(int i = 0; i < rows; i++){

            // Validate A[i][i] is not zero

            if( A->data[ i*(rows+1) ] == 0 ){

                printf("\nDivision by 0. Cannot continue. Stopped at iteration: %d\n\n", iter+1);
            }

            sum = 0;

            for(int j = 0; j < rows; j++){

                if(j != i)
                    sum += A->data[ i*rows + j ] * x_new->data[j];
            }

            temp_x = ( b->data[i] - sum ) / A->data[ i*(rows+1) ];

            num_norm += (temp_x - x_new->data[i]) * (temp_x - x_new->data[i]);
            x_new->data[i] = temp_x;
            den_norm += x_new->data[i] * x_new->data[i];
        }

        error = sqrt(num_norm/den_norm);

        // If convergence is reached, stop iterating and return solution x_new

        if( error < tolerance ){
            printf("\n***********\n\nConvergence reached at iteration %d\n\n***********\n", iter+1);
            return x_new;
        }

    }

    printf("\nSolution x did not converge\n");

    return x_new;
}

double bisectionSolver( double (*f)(double), double xmin, double xmax, double tol ){

    if( f(xmin)*f(xmax) >= 0 ){
        printf("\nBisection method fails\n");
        exit(-1);
    }

    double xmid;

    while ( (xmax - xmin) > tol ) {

        xmid = (xmin + xmax) / 2;

        if ( f(xmid) * f(xmin) < 0 ) {

            xmax = xmid;

        } else if ( f(xmid) * f(xmax) < 0 ){

            xmin = xmid;

        } else if ( f(xmid) * f(xmax) < 0 ){

            return xmid;

        } else {

            printf("\nBisection method fails\n");
            exit(-1);

        }
    }

    return (xmin + xmax) / 2;
}

double newtonSolver( double (*f)(double), double (*fp)(double), double x0, double tol, double epsilon, int numIter ){

    double xn;

    for (int i = 0; i < numIter; i++) {

        double y  = f(x0);
        double yp  = fp(x0);

        if ( fabs( yp ) <= epsilon ) {
            printf("\nSolution did not converge\n");
            exit(-1);
        }

        xn = x0 - y/yp;

        if ( fabs( xn - x0 ) < tol ) {

            printf("\nSolution converged at iteration: %d\n", i+1);
            return xn;
        }

        x0 = xn;
    }

    printf("Solution did not converge\n");

    return xn;
}

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

double * powerSolver(Matrix * A, Matrix * eigenVectOld, double* lambdaInit, int num_iters, double epsilon){

    // Size of matrix and vectors
    int rows = A->rows;

    // Variable to update dominant eigenvalue
    double lambdaNew;

    // Helper variable for dot product of vectors
    double accum;

    // Array to plot convergence of eigenvalue calculated
    double* lambdaTrace = malloc( num_iters * sizeof( *lambdaTrace ) );
    int lambdaTraceSize = 0;

    // Variable to iterate algorithm
    int i = 0;

    // Vector to calculate A * v_(k-1)
    Matrix *eigenVectNew = malloc( sizeof( eigenVectNew ) );

    // Normalize initial vector v_0
    double norm = vectNorm(eigenVectOld);
    for (int k = 0; k < rows; k++)
        eigenVectOld->data[k] *= (1/norm);

    // Start iterations

    for (i = 0; i < num_iters; i++) {

        // Calculate w =  A * v_(k-1)
        eigenVectNew = multiply(A, eigenVectOld);

        // Calculate v_(K) = w/norm(w)
        norm = vectNorm(eigenVectNew);
        for (int k = 0; k < rows; k++)
            eigenVectOld->data[k] = eigenVectNew->data[k] * (1/norm);

        // Calculate dominant eigenvalue and store in lambdaNew

        accum = 0;

        for (int j = 0; j < rows; j++) {
            accum += eigenVectNew->data[j] * eigenVectOld->data[j];
        }

        lambdaNew = accum;

        // Check for convergence and stop or update the eigenvalue

        if ( fabs( lambdaNew - (*lambdaInit) ) / lambdaNew < epsilon ) {
            (*lambdaInit) = lambdaNew;
            break;
        }
        else{
            lambdaTrace[i] = lambdaNew;
            (*lambdaInit) = lambdaNew;
        }
    }

    printf("\nConverged after %d iterations\n\n", i);

    return lambdaTrace;
}
