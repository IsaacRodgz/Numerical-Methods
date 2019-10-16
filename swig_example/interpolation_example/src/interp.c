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

void newton_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis){

    double* coeffs = malloc( n * sizeof *coeffs );

    for (int i = 0; i < n; i++) {
        coeffs[i] = fis[i];
    }

    for (int i = 1; i < n; i++) {

        for (int j = 0; j <= i-1 ; j++) {

            coeffs[i] = ( coeffs[j] - coeffs[i] ) / ( xis[j] - xis[i] );
        }
    }

    for (int i = 0; i < m; i++) {

        yis[i] = 0.0;
        double x_temp = 1.0;

        for(int j = 0; j < n; j++){

            yis[i] += coeffs[j]*x_temp;
            x_temp *= (pts[i] - xis[j]);
        }
    }

    free(coeffs);
}
