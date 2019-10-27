/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-10-15
*/

#ifndef _integration_h
#define _integration_h
#include "integration.h"

double trapezoidal_rule( double (* f)(double), double x_left, double x_right, int interval_size );

double simpson3_rule( double (* f)(double), double x_left, double x_right, int interval_size );

double simpson8_rule( double (* f)(double), double x_left, double x_right, int interval_size );

double boole_rule( double (* f)(double), double x_left, double x_right, int interval_size );

double boole_rule( double (* f)(double), double x_left, double x_right, int interval_size );

double weddle_rule( double (* f)(double), double x_left, double x_right, int interval_size );

/*
Test mathematical functions for integration
*/

double sin_x(double x);

double sin_sqr_x(double x);

#endif
