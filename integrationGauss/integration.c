/*
Author: Isaac Rodríguez Bribiesca
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "integration.h"
#define TRUE 1
#define FALSE 0

double one_point_gl( double (* f)(double) ){

    return 2*f(0);
}

double two_point_gl( double (* f)(double) ){

    return (f(-1/sqrt(3)) + f(1/sqrt(3)));
}

double three_point_gl( double (* f)(double) ){

    return ( (5.0/9.0)*f(-sqrt(3.0/5.0)) + (8.0/9.0)*f(0) + (5.0/9.0)*f(sqrt(3.0/5.0)) );
}

double four_point_gl( double (* f)(double) ){

    double a0 = (18.0+sqrt(30))/36.0;
    double a1 = (18.0-sqrt(30))/36.0;
    double x0 = sqrt((3-2*sqrt(6.0/5.0))/7.0);
    double x1 = sqrt((3+2*sqrt(6.0/5.0))/7.0);

    return (a0*f(x0) + a0*f(-x0) + a1*f(x1) + a1*f(-x1));
}

/*
Test mathematical functions for integration
*/

double f1(double x) {

    return 1.0/(1.0+(x*x));
}

double f2(double x) {

    return 0.5*sqrt( 1 + cos( 0.5*x + 1.5 )*cos( 0.5*x + 1.5 ) );
}
