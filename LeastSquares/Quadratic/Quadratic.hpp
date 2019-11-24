#ifndef QUADRATIC_HPP
#define QUADRATIC_HPP

#include <bits/stdc++.h>
#include <string.h>

using namespace std;

class Quadratic {

public:

    vector<double> coeffs;
    vector<double> b;
    vector<vector<double> > A;

    vector<double> x;
    vector<double> y;
    int size;

    Quadratic(int sizep);

    double f_quadratic(double x);

    void create_quadratic_data();

    void fill_quadratic_system();

    void solve_quadratic();

};

#endif
