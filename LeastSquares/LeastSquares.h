#ifndef LEAST_SQUARES_H
#define LEAST_SQUARES_H

#include <bits/stdc++.h>
#include <string.h>

using namespace std;

class LeastSquares {

public:

    vector<double> coeffs;
    vector<double> b;
    vector<vector<double> > A;

    vector<double> x;
    vector<double> y;
    int size;

    LeastSquares(int sizep);

    double f(double x);

    void create_data();

    void fill_system();

    void solve();

};

#endif
