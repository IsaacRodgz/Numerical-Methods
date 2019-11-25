#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP

#include <bits/stdc++.h>
#include <string.h>

using namespace std;

class FiniteElement {

public:

    vector<double> nodes;
    vector<double> coeffs;
    vector<double> b;
    vector<vector<double> > A;

    vector<double> x;
    vector<double> y;

    double lambda;
    int num_elements;
    int size;

    FiniteElement(double lambda_p, int num_elements_p);

    double f_linear(double x);

    void create_linear_data();

    void read_linear_data(string data_file);

    double map(double t, double xmin, double xmax);

    double n_i1(double r);

    double n_i2(double r);

    void fill_fe_system();

    void solve_fem();

};

#endif
