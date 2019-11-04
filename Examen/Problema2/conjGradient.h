#ifndef conjGradient_H
#define conjGradient_H

#include<bits/stdc++.h>

using namespace std;

int is_simetric(vector< vector<double> > A);

vector<double> read_vector(const string vector_file);

vector< vector<double> > read_matrix(const string matrix_file);

vector< vector<double> > read_matrix_full(const string matrix_file);

tuple< vector<double>, vector<int>, vector<int> > read_matrix_rala(const string matrix_file);

vector<double> cGradientSolver(const string matrix_file, const string vector_file, int numIters, double epsilon);

#endif
