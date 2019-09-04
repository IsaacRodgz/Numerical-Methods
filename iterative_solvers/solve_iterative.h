/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-26
Description: This program includes declaration of functions that solve systems A*x = b of n equations with n variables,
through two different iterative methodos methods, namely, Jacobi and Gauss-Seidel.
*/

#ifndef _solve_iterative_h
#define _solve_iterative_h

#include "matrix_struct.h"

Matrix * jacobiSolver(Matrix * A, Matrix * b, int num_iters, double tolerance);

Matrix * gaussSeidelSolver(Matrix * A, Matrix * b, int num_iters, double tolerance);

double vectNorm(Matrix* vect);

void swap(Matrix ** A, Matrix ** B);

double * powerSolver(Matrix * A, Matrix * eigenVectOld, double* lambdaInit, int num_iters, double epsilon);

double bisectionSolver( double (*f)(double), double xmin, double xmax, double tol, int numIters );

double newtonSolver( double (*f)(double), double (*fp)(double), double x0, double tol, double epsilon, int numIter );

#endif
