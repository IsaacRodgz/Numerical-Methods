/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-10-15
*/

#ifndef _interp_h
#define _interp_h
#include "interp.h"

void newton_interp(int n, double* xis, int m, double* fis, int o, double* coeffs);

double newton_eval(int n, double* coeffs, int m, double* xis, double x);

#endif
