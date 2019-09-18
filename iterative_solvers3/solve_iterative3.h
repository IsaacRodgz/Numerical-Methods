/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-26
Description: This program includes declaration of functions that solve systems A*x = b of n equations with n variables,
through two different iterative methodos methods, namely, Jacobi and Gauss-Seidel.
*/

#ifndef _solve_iterative3_h
#define _solve_iterative3_h

#include "matrix_struct.h"

void rayleighSolver(Matrix * A, Matrix * eigenVect, double * lambda, int num_iters, double epsilon);

void subspaceSolver(Matrix * A, Matrix * FI, Matrix * LA, int num_iters, double epsilon, int k);

void cGradientSolver(Matrix * A, Matrix * b, Matrix * x, int numIters, double epsilon);

void QRFactor(Matrix * A, Matrix * Q, Matrix * R, int numIters, double epsilon);

#endif
