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

int is_simetric(Matrix *A);

void kPowerSolver(Matrix * A, Matrix * eigenVects, Matrix * eigenVals, int num_iters, double epsilon, int k);

void kInversePowerSolver(Matrix * A, Matrix * eigenVects, Matrix * eigenVals, int num_iters, double epsilon, int k);

int is_diagonal(Matrix *A);

Matrix * solve_lower_triang(Matrix *A, Matrix *b, int fill_diag);

Matrix * solve_upper_triang(Matrix *A, Matrix *b, int fill_diag);

void set_pivot_row(Matrix *A, int* pivot, int i_switch, int index_set);

int max_pivot_index_column(Matrix *A, int limit);

void factor_doolittle_pivoting(Matrix *A, int *pivot);

Matrix * solve_doolittle_pivoting(Matrix *A, Matrix *b, int *pivots);

#endif
