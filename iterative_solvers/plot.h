/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-26
Description: This program includes declaration of functions that solve systems A*x = b of n equations with n variables,
through two different iterative methodos methods, namely, Jacobi and Gauss-Seidel.
*/

#ifndef _plot_h
#define _plot_h

#include "plot.h"

void plotFunction( double (*f)(double), double xmin, double xmax, int nIntervals, double xSolution);

void plotData( double * data, int size );

double f1(double x);

double f2(double x);

double f3(double x);

double f4(double x);

double f5(double x);

double f6(double x);

#endif
