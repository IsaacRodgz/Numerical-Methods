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

#endif
