/*
Author: Isaac Rodríguez Bribiesca
*/

#ifndef _integration_h
#define _integration_h
#include "integration.h"

double one_point_gl( double (* f)(double) );

double two_point_gl( double (* f)(double) );

double three_point_gl( double (* f)(double) );

double four_point_gl( double (* f)(double) );

/*
Test mathematical functions for integration
*/

double f1(double x);

double f2(double x);

#endif
