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

        double s = (pts[i] - xis[0])/(xis[1]-xis[0]);
        double s_accum = s;
        double fact = 1.0;
        yis[i] = coeffs[0] + (s*coeffs[1])/fact;

        for(int j = 2; j < n; j++){

            s_accum *= (s-j+1);
            fact *= j;
            yis[i] += (s_accum*coeffs[j])/fact;
        }
    }

    free(coeffs);
}

void gregory_backward_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis) {

    double* coeffs = malloc( n * sizeof *coeffs );

    for (int i = 0; i < n; i++) {
        coeffs[i] = fis[i];
    }

    for (int i = 1; i < n; i++) {

        for (int j = 0; j < n-i ; j++) {

            coeffs[j] = coeffs[j+1] - coeffs[j];
        }
    }

    for (int i = 0; i < m; i++) {

        double s = (pts[i] - xis[n-1])/(xis[1]-xis[0]);
        double s_accum = s;
        double fact = 1.0;
        yis[i] = coeffs[n-1] + (s*coeffs[n-2])/fact;

        for(int j = 2; j < n; j++){

            s_accum *= (s+j-1);
            fact *= j;
            yis[i] += (s_accum*coeffs[n-j-1])/fact;
        }
    }

    free(coeffs);
}
