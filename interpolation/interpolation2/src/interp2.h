/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-10-15
*/

#ifndef _interp2_h
#define _interp2_h
#include "interp2.h"

void gregory_forward_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis);

void gregory_backward_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis);

/*
void newton_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis);

void lagrange_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis);

void hermite_interp_t(int n, double* xis, double* fis,  double* fpis, int m, double* pts, double* yis);

void newton_piecewise_interp_t(int n, double* xis, double* fis, int m, double* pts, double* yis);
*/

#endif
