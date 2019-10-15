/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-10-15
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "interp.h"
#define TRUE 1
#define FALSE 0

void newton_interp(double* xis, double* fis, int n ){

    double* coeffs = malloc( n * sizeof *coeffs );

    for (int i = 0; i < n; i++) {
        coeffs[i] = fis[i];
    }

    for (int i = 1; i < n; i++) {

        for (int j = 0; j <= i-1 ; j++) {

            coeffs[i] = ( coeffs[j] - coeffs[i] ) / ( xis[j] - xis[i] );

        }
    }

    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("%lf ", coeffs[i]);
    }
    printf("\n\n");
}

double newton_eval(){}
