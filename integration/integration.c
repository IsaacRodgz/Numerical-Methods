/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-10-15
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "integration.h"
#define TRUE 1
#define FALSE 0

double trapezoidal_rule( double (* f)(double), double x_left, double x_right, int interval_size ) {

    double delta = (x_right - x_left)/interval_size;

    double integral = f(x_left) +  f(x_right);

    for (int i = 1; i < interval_size; i++) {

        integral += 2*f(x_left + i*delta);
    }

    integral *= (delta/2);

    return integral;
}

double simpson3_rule( double (* f)(double), double x_left, double x_right, int interval_size ) {

    double delta = (x_right - x_left)/interval_size;

    double integral = f(x_left) +  f(x_right);

    for (int i = 1; i < interval_size; i++) {

        if ( i%2 == 0 )
            integral += 2*f(x_left + i*delta);
        else
            integral += 4*f(x_left + i*delta);

    }

    integral *= (delta/3);

    return integral;
}

double simpson8_rule( double (* f)(double), double x_left, double x_right, int interval_size ) {

    double delta = (x_right - x_left)/interval_size;

    double integral = f(x_left) +  f(x_right);

    //integral += 3*f(x_left + delta);
    //integral += 3*f(x_left + 2*delta);

    for (int i = 1; i < interval_size; i++) {

        if ( i%3 == 0 )
            integral += 2*f(x_left + i*delta);
        else
            integral += 3*f(x_left + i*delta);
    }

    integral *= (delta*(3.0/8.0));

    return integral;
}

double boole_rule( double (* f)(double), double x_left, double x_right, int interval_size ) {

    double delta = (x_right - x_left)/(4*interval_size);

    double integral = 7*f(x_left);

    double x_curr = x_left + delta;

    for (int i = 1; i < interval_size; i++) {

        integral += 32*f(x_curr);
        x_curr += delta;
        integral += 12*f(x_curr);
        x_curr += delta;
        integral += 32*f(x_curr);
        x_curr += delta;
        integral += 14*f(x_curr);
        x_curr += delta;
    }

    integral += 32*f(x_curr);
    x_curr += delta;
    integral += 12*f(x_curr);
    x_curr += delta;
    integral += 32*f(x_curr);
    integral += 7*f(x_right);

    integral *= (delta*(2.0/45.0));

    return integral;
}

double weddle_rule( double (* f)(double), double x_left, double x_right, int interval_size ) {

    double delta = (x_right - x_left)/(6*interval_size);

    double integral = 41*f(x_left);

    double x_curr = x_left + delta;

    for (int i = 1; i < interval_size; i++) {

        integral += 216*f(x_curr);
        x_curr += delta;
        integral += 27*f(x_curr);
        x_curr += delta;
        integral += 272*f(x_curr);
        x_curr += delta;
        integral += 27*f(x_curr);
        x_curr += delta;
        integral += 216*f(x_curr);
        x_curr += delta;
        integral += 82*f(x_curr);
        x_curr += delta;
    }

    integral += 216*f(x_curr);
    x_curr += delta;
    integral += 27*f(x_curr);
    x_curr += delta;
    integral += 272*f(x_curr);
    x_curr += delta;
    integral += 27*f(x_curr);
    x_curr += delta;
    integral += 216*f(x_curr);
    integral += 41*f(x_right);

    integral *= (delta/140.0);

    return integral;
}

/*
Test mathematical functions for integration
*/

double sin_x(double x) {

    return sin(x + 0.5);
}

double sin_sqr_x(double x) {

    return sin(x)*sin(x);
}

double f1(double x) {

    return 1.0/(1.0+(x*x));
}

double f2(double x) {

    return sqrt( 1 + cos(x)*cos(x) );
}
