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

void newton_interp(int n, double* xis, int m, double* fis, int o, double* coeffs){

    //double* coeffs = malloc( n * sizeof *coeffs );

    for (int i = 0; i < n; i++) {
        coeffs[i] = fis[i];
    }

    for (int i = 1; i < n; i++) {

        for (int j = 0; j <= i-1 ; j++) {

            coeffs[i] = ( coeffs[j] - coeffs[i] ) / ( xis[j] - xis[i] );
        }
    }
}

double newton_eval(int n, double* coeffs, int m, double* xis, double x){

    double pn = 0.0;
    double x_temp = 1.0;

    for(int i = 0; i < n; i++){
        
        pn += coeffs[i]*x_temp;
        x_temp *= (x - xis[i]);
    }

    return pn;
}

