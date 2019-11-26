#include <math.h>
#include <bits/stdc++.h>

#include "Spectral.hpp"

using namespace std;

Spectral::Spectral(double sigma_p) : sigma(sigma_p) {}

void Spectral::readData(int num_cols, string data_file){

    ifstream file(data_file);
    string line;

    while ( getline(file, line) ) {

        stringstream ss(line);

        vector<double> point;
        double elem;

        for (int i = 0; i < num_cols; i++) {

            ss >> elem;

            point.push_back(elem);
        }

        data.push_back(point);
    }

    file.close();
}

double Spectral::gaussianDistance(vector<double> x1, vector<double> x2){

    vector<double> x3;

    for (int i = 0; i < x1.size(); i++) {

        double diff = x1[i] - x2[i];

        x3.push_back(diff);
    }

    return exp(-norm_squared(x3)/(2*sigma));
}

double Spectral::norm_squared(vector<double> x){

    double norm_val = 0;

    for (int i = 0; i < x.size(); i++) {

        norm_val += x[i]*x[i];
    }

    return norm_val;
}

void Spectral::buildGraph(){

    graph.resize(data.size(), vector<double>(data.size()));

    for (int i = 0; i < graph.size(); i++) {

        for (int j = 0; j < graph.size(); j++) {

            graph[i][j] = gaussianDistance(data[i], data[j]);
        }
    }
}

void Spectral::buildLaplacian(){

    vector<double> diag;

    for (int i = 0; i < graph.size(); i++) {

        double sum = 0;

        for (int j = 0; j < graph.size(); j++) {

            sum += graph[i][i];
        }

        diag.push_back(sum);
    }

    lapacian.resize(graph.size(), vector<double>(graph.size()));

    for (int i = 0; i < graph.size(); i++) {

        for (int j = 0; j < graph.size(); j++) {

            if ( i == j ) {

                lapacian[i][j] = graph[i][j] - diag[i];
            }

            else {

                lapacian[i][j] = graph[i][j];
            }
        }
    }
}

void Spectral::computeEigen(){

    
}

// Print functions

void Spectral::printData(){

    for (int i = 0; i < data.size(); i++) {

        for (int j = 0; j < data[0].size(); j++) {

            cout << data[i][j] << ", ";
        }

        cout << endl;
    }
}

void Spectral::printGraph(){

    cout << endl;

    for (int i = 0; i < graph.size(); i++) {

        for (int j = 0; j < graph[0].size(); j++) {

            printf("%10lf ", graph[i][j]);
        }

        cout << endl;
    }
}
