/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-26
Description: This program includes declaration of functions that solve systems A*x = b of n equations with n variables,
through two different iterative methodos methods, namely, Jacobi and Gauss-Seidel.
*/

#ifndef _solve_iterative_h
#define _solve_iterative_h

#include "matrix_struct.h"

double vectNorm(Matrix* vect);

void swap(Matrix ** A, Matrix ** B);

void powerSolver(Matrix * A, Matrix * eigenVectOld, double* lambdaInit, int num_iters, double epsilon);

void kPowerSolver(Matrix * A, Matrix * eigenVects, Matrix * eigenVals, int num_iters, double epsilon, int k);

void inversePowerSolver(Matrix * A, Matrix * eigenVect, double* lambdaInit, int num_iters, double epsilon);

void kInversePowerSolver(Matrix * A, Matrix * eigenVects, Matrix * eigenVals, int num_iters, double epsilon, int k);

void maxOffDiagonal(Matrix * A, int* p, int* q);

void jacobiSolver(Matrix * A, Matrix * F, int num_iters, double epsilon);

#endif