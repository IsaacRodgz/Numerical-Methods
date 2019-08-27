/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-26
Description: This program includes declaration of functions that solve systems A*x = b of n equations with n variables,
through two different iterative methodos methods, namely, Jacobi and Gauss-Seidel.
*/

#ifndef _plot_h
#define _plot_h

#include "plot.h"

void plotData( double (*f)(double), double xmin, double xmax, int nIntervals, double xSolution);

double polynomial4(double x);

double polynomial4P(double x);

double squareRoot(double x);

double squareRootP(double x);

#endif
