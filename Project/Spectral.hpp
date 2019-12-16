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

    vector< vector<double> > eigenData;

    vector<int> clusterIDs;

    vector<vector<double> > centers;

    vector< vector<double> > graph;

    vector< vector<double> > laplacian;

    vector< vector<double> > eigenVects;

    vector< vector<double> > eigenVals;
    vector<double> eigenVals_vect;

    // Cairo graphing method

    void plot();

    // Graph and laplacian construction methods

    void readData(int num_cols, string data_file);

    double gaussianDistance(vector<double> x1, vector<double> x2);

    double norm_squared(vector<double> x);

    void buildGraph();

    void buildLaplacian();

    // Print functions

    void printData();

    void printGraph();

    void printClusters();

    void printEigen();

    void printLaplacian();

    // k-means methods

    void kMeans(int numIters, int numClusters);

    void kMeansEig(int numIters, int numClusters);

    int getNearestCenter(vector<double> p);

    // Eigensolver methods

    int is_diagonal(vector< vector<double> > A, double epsilon);

    void computeEigen(int numEigen, int numIters, double epsilon);

    void transformData(int numEigen);

    void subspaceSolver(int num_iters, double epsilon, int k);

    void QRFactor(vector< vector<double> > &A, vector< vector<double> > &Q, vector< vector<double> > &R);

    vector<double> solve_lower_triang(vector< vector<double> > &A, vector<double> &b, int fill_diag);

    vector<double> solve_upper_triang(vector< vector<double> > &A, vector<double> &b, int fill_diag);

    void set_pivot_row(vector< vector<double> > &A, vector<int> &pivot, int i_switch, int index_set);

    int max_pivot_index_column(vector< vector<double> > &A, int limit);

    void factor_doolittle_pivoting(vector< vector<double> > &A, vector<int> &pivot);

    void factor_doolittle(vector< vector<double> > &A);

    vector<double> solve_doolittle_pivoting(vector< vector<double> > &A, vector<double> &b, vector<int> &pivots);

    vector<double> solve_doolittle(vector< vector<double> > &A, vector<double> &b, int factor_flag);

    vector< vector<double> > solve_inverse(vector< vector<double> > &A);

    int is_simetric(vector< vector<double> > &A);

    void deflation(vector< vector<double> > &eigenVects, vector<double> &eigenVectInit, int currCol);

    void kInversePowerSolver(vector< vector<double> > &A, vector< vector<double> > &eigenVects, vector<double> &eigenVals, int num_iters, double epsilon, int k);

    double vectNorm(vector<double> &vect);

    void swap(vector<double> &A, vector<double> &B);

    void maxOffDiagonal(vector< vector<double> > &A, int &p, int &q);

    void jacobiSolver(vector< vector<double> > &A, vector< vector<double> > &F, int num_iters, double epsilon);
};

#endif
