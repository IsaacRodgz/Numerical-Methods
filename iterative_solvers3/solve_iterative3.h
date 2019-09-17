/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-26
Description: This program includes declaration of functions that solve systems A*x = b of n equations with n variables,
through two different iterative methodos methods, namely, Jacobi and Gauss-Seidel.
*/

#ifndef _solve_iterative3_h
#define _solve_iterative3_h

#include "matrix_struct.h"

double vectNorm(Matrix* vect);

void rayleighSolver(Matrix * A, Matrix * eigenVect, double * lambda, int num_iters, double epsilon);

#endif