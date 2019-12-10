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

    vector< vector<double> > laplacian;

    vector< vector<double> > eigenVects;

    vector< vector<double> > eigenVals;

    void readData(int num_cols, string data_file);

    double gaussianDistance(vector<double> x1, vector<double> x2);

    double norm_squared(vector<double> x);

    void buildGraph();

    // Print functions

    void printData();

    void printGraph();

    void buildLaplacian();

    int is_diagonal(vector< vector<double> > A);

    void computeEigen(int numEigen, int numIters, double epsilon);

    void subspaceSolver(int num_iters, double epsilon, int k);

    void QRFactor(vector< vector<double> > &A, vector< vector<double> > &Q, vector< vector<double> > &R);

};

#endif
