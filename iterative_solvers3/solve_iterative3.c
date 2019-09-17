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
