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

            num_norm += pow( x_new->data[i] - x_old->data[i] , 2);
            den_norm += pow(x_new->data[i], 2);
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

                if(j < i)
                    sum += A->data[ i*rows + j ] * x_new->data[j];
                if(j > i)
                    sum += A->data[ i*rows + j ] * x_old->data[j];
            }

            x_new->data[i] = ( b->data[i] - sum ) / A->data[ i*(rows+1) ];

            num_norm += pow( x_new->data[i] - x_old->data[i] , 2);
            den_norm += pow(x_new->data[i], 2);
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

void normalizeVector(Matrix * vect){

    double sum = 0;

    for (int i = 0; i < vect->rows; i++) {
        sum += vect->data[i];
    }

    if ( (sum = sqrt(sum)) == 0.0 ) {
        fprintf(stderr, "\n[Error] Vector of norm 0 found\n\n");
    }

    for (int i = 0; i < vect->rows; i++) {
        vect->data[i] *= (1/sum);
    }
}

double * powerSolver(Matrix * A, Matrix * eigenVectOld, double* lambdaInit, int num_iters, double tolerance, double epsilon){

    int rows = A->rows;
    double lambdaNew;
    double numerator;
    double denominator;
    double* lambdaTrace = malloc( num_iters * sizeof( *lambdaTrace ) );
    int lambdaTraceSize = 0;

    Matrix *eigenVectNew = malloc( sizeof( eigenVectNew ) );

    for (int i = 0; i < num_iters; i++) {

        normalizeVector(eigenVectOld);
        eigenVectNew = multiply(A, eigenVectOld);

        numerator = 0;
        denominator = 0;

        for (int j = 0; j < rows; j++) {

            numerator += eigenVectNew->data[j] * eigenVectNew->data[j];
            denominator += eigenVectNew->data[j] * eigenVectOld->data[j];
        }

        if ( denominator > tolerance )
            lambdaNew = numerator / denominator;
        else{
            fprintf(stderr, "\n[Error]: denominator close or equal to zero\n\n", );
        }

        if ( fabs( lambdaNew - lambdaInit ) < epsilon ) {
            *lambdaInit = lambdaNew;
            break;
        }
        else{
            lambdaTrace[i] = lambdaNew
            *lambdaInit = lambdaNew;
            // TODO: SWAP( eigenVectNew, eigenVectOld )
        }
    }
    // TODO: Plot lambda values
    return lambdaTrace;
}
