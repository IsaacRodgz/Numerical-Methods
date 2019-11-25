#ifndef LINEAR_HPP
#define LINEAR_HPP

#include <bits/stdc++.h>
#include <string.h>

using namespace std;

class Linear {

public:

    vector<double> coeffs;
    vector<double> b;
    vector<vector<double> > A;

    vector<double> x;
    vector<double> y;
    int size;

    double f_linear(double x);

    void create_linear_data();

    void read_linear_data(string data_file);

    void fill_linear_system();

    void solve_linear();

};

#endif
