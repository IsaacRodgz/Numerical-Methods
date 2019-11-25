#include "Linear.hpp"

#include <bits/stdc++.h>
#include <string.h>
#include <random>

using namespace std;

// Linear model

double Linear::f_linear(double x){

    return 1.37*x + 2.43;
}

void Linear::create_linear_data(){

    size = 32;

    default_random_engine generator;
    random_device rd;
    mt19937 gen(rd());

    double xi = -5;
    double step = 10.0/double(size);
    int i = 0;

    while (xi < 5) {

        x.push_back(xi);
        normal_distribution<double> distribution( f_linear(x[i]) , 1.0 );
        y.push_back(distribution(generator));

        xi += step;
        i++;
    }
}

void Linear::read_linear_data(string data_file){

    ifstream file(data_file);
    string line;

    double val_x;
    double val_y;

    while ( getline(file, line) ) {

        stringstream ss(line);

        ss >> val_x >> val_y;

        x.push_back(val_x);
        y.push_back(val_y);

    }

    file.close();
    size = x.size();
}

void Linear::fill_linear_system(){

    A.resize(2, vector<double>(2, 0));
    b.resize(2, 0);

    for (int i = 0; i < size; i++) {

        A[0][0] += x[i]*x[i];
        A[0][1] += x[i];
        A[1][0] += x[i];

        b[0] += x[i]*y[i];
        b[1] += y[i];
    }

    A[1][1] = size;
}

void Linear::solve_linear(){

    coeffs.resize(2, 0);
    vector<vector<double> > A_inv(2, vector<double>(2, 0));
    double det = 0;

    det = 1.0/(A[0][0]*A[1][1] - A[0][1]*A[1][0]);

    A_inv[0][0] = det*A[1][1];
    A_inv[0][1] = -det*A[0][1];
    A_inv[1][0] = -det*A[1][0];
    A_inv[1][1] = det*A[0][0];

    for (int i = 0; i < 2; i++) {

        double sum = 0;

        for (int j = 0; j < 2; j++) {

            sum += A_inv[i][j]*b[j];
        }

        coeffs[i] = sum;
    }
}
