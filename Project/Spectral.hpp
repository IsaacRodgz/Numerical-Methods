#ifndef SPECTRAL_HPP
#define SPECTRAL_HPP

#include <bits/stdc++.h>
#include <string.h>

using namespace std;

class Spectral {

public:

    Spectral(double sigma_p);

    double sigma;

    vector< vector<double> > data;

    vector< vector<double> > graph;

    vector< vector<double> > lapacian;

    void readData(int num_cols, string data_file);

    double gaussianDistance(vector<double> x1, vector<double> x2);

    double norm_squared(vector<double> x);

    void buildGraph();

    // Print functions

    void printData();

    void printGraph();

    void buildLaplacian();

    void computeEigen();

    void subspaceSolver(vector< vector<double> > A, vector< vector<double> > FI, vector< vector<double> > LA, int num_iters, double epsilon, int k);

    void QRFactor(vector< vector<double> > A, vector< vector<double> > Q, vector< vector<double> > R);

};

#endif
