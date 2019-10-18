/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-10-15
*/

#ifndef _interp_h
#define _interp_h
#include "interp.h"

void newton_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis);

void lagrange_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis);

void hermite_interp_t(int n, double* xis, double* fis,  double* fpis, int m, double* pts, double* yis);

#endif
