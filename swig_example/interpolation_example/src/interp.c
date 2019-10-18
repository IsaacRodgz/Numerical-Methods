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

void lagrange_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis){

    for (int k = 0; k < m; k++) {

        yis[k] = 0;

        for (int i = 0; i < n; i++) {

            double term = fis[i];

            for (int j = 0; j <= n ; j++) {

                if ( i != j ) {
                    term = term*(pts[k] - xis[j])/(xis[i] - xis[j]);
                }
            }

            yis[k] += term;
        }
    }
}

void hermite_interp_t(int n, double* xis, double* fis,  double* fpis, int m, double* pts, double* yis){

    double** diffs = malloc( 2*n * sizeof *diffs );
    diffs[0] = malloc( 4*n*n * sizeof **diffs );
    for (int i = 1; i < 2*n; i++) {
        diffs[i] = diffs[i-1] + 2*n;
    }

    double* z = malloc( 2*n * sizeof *diffs );

    // Initialize z_(2i) = x_i and z_(2i+1) = x_i
    // And f[z_(2i)] = f(x_i), f[z_(2i+1)] = f(x_i)

    for (int i = 0; i < n; i++) {

        z[2*i] = xis[i];
        z[2*i+1] = xis[i];

        diffs[2*i][0] = fis[i];
        diffs[2*i+1][0] = fis[i];

        diffs[2*i+1][1] = fpis[i];

        if ( i != 0 )
            diffs[2*i][1] = (diffs[2*i][0] - diffs[2*i-1][0]) / (z[2*i] - z[2*i-1]);
    }

    for (int i = 2; i < 2*n; i++) {

        for (int j = 2; j <= i ; j++) {

            diffs[i][j] = (diffs[i][j-1] - diffs[i-1][j-1]) / (z[i] - z[i-j]);
        }
    }

    printf("\n");
    for (int i = 0; i < 2*n; i++) {
        for (int j = 0; j < 2*n; j++){
            printf("%10lf", diffs[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    for (int i = 0; i < m; i++) {

        yis[i] = diffs[0][0];
        double x_temp = 1.0;

        for(int j = 1; j < n; j++){

            x_temp *= (pts[i] - xis[j-1]);
            yis[i] += diffs[2*j-1][2*j-1]*x_temp;

            x_temp *= (pts[i] - xis[j-1]);
            yis[i] += diffs[2*j][2*j]*x_temp;

        }
    }

    free(diffs[0]);
    free(diffs);
    free(z);
}
