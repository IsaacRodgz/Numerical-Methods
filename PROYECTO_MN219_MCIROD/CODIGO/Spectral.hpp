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

    void kMeansEig(int numIters, int numClusters);

    int getNearestCenter(vector<double> p);

    // Eigensolver methods

    int is_diagonal(vector< vector<double> > A, double epsilon);

    void computeEigen(int numEigen, int numIters, double epsilon);

    void transformData(int numEigen);

    int is_simetric(vector< vector<double> > &A);

    double vectNorm(vector<double> &vect);

    void maxOffDiagonal(vector< vector<double> > &A, int &p, int &q);

    void jacobiSolver(vector< vector<double> > &A, vector< vector<double> > &F, int num_iters, double epsilon);
};

#endif
