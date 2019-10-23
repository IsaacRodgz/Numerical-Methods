/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-10-15
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "interp2.h"
#define TRUE 1
#define FALSE 0

void gregory_forward_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis) {

    double* coeffs = malloc( n * sizeof *coeffs );

    for (int i = 0; i < n; i++) {
        coeffs[i] = fis[i];
    }

    for (int i = 1; i < n; i++) {

        for (int j = n-1; j >= i ; j--) {

            coeffs[j] -= coeffs[j-1];
        }
    }

    for (int i = 0; i < m; i++) {

        double s = (pts[i] - xis[0])/(fabs(xis[0]-xis[1]));
        double bin_coeff = 

        yis[i] = 0.0;


        for(int j = 0; j < n; j++){

            yis[i] += coeffs[j]*x_temp;
            x_temp *= (pts[i] - xis[j]);
        }
    }

    free(coeffs);
}
