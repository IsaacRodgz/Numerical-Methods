/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-12
*/

#ifndef _solve_matrix_direct_h
#define _solve_matrix_direct_h

#include "matrix_struct.h"

int is_diagonal(Matrix *A);

Matrix * solve_lower_triang(Matrix *A, Matrix *b, int fill_diag);

Matrix * solve_upper_triang(Matrix *A, Matrix *b, int fill_diag);

void set_pivot_row(Matrix *A, int* pivot, int i_switch, int index_set);

int max_pivot_index_column(Matrix *A, int limit);

void factor_doolittle(Matrix *A);

void factor_doolittle_pivoting(Matrix *A, int *pivot);

Matrix * solve_doolittle(Matrix *A, Matrix *b, int factor_flag);

Matrix * solve_doolittle_pivoting(Matrix *A, Matrix *b, int *pivots);

Matrix * solve_inverse(Matrix *A);

#endif
