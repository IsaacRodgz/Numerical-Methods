/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-12
Description: Header File for pogram that solves system A*x = b of n equations with n variables
through different methods.

- Diagonal matrix - solve_diagonal
- Lower triangular matrix - solve_lower_triang
- Upper triangular - solve_upper_triang
- Gauss Elimination, assume pivot != 0. Reduces A to upper triangular - solve_gauss_elim

*/

#ifndef _solve_matrix_direct_h
#define _solve_matrix_direct_h

#include "matrix_struct.h"

int test_diagonal(Matrix *A);

Matrix * solve_diagonal(Matrix *A, Matrix *b);

Matrix * solve_lower_triang(Matrix *A, Matrix *b, int fill_diag);

Matrix * solve_upper_triang(Matrix *A, Matrix *b, int fill_diag);

double diagonal_determinant(Matrix *A);

void set_pivot(Matrix *A, Matrix *b, int *index_order, int i_switch, int j_switch, int index_set);

int * max_pivot_index(Matrix *A, int limit, int *pivot_index);

void solve_gauss_elim(Matrix *A, Matrix *b);

void solve_gauss_elim_pivot(Matrix *A, Matrix *b, int *index_order);

Matrix * solve_no_pivot(Matrix *A, Matrix *b);

Matrix * solve_pivot(Matrix *A, Matrix *b);

void factor_doolittle(Matrix *A);

Matrix * solve_doolittle(Matrix *A, Matrix *b, int factor_flag);

Matrix * solve_inverse(Matrix *A);

void factor_cholesky(Matrix *A);

Matrix * solve_cholesky_modified(Matrix *A, Matrix *b);

#endif
