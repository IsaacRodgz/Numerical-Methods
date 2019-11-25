#include "Quadratic.hpp"

#include <bits/stdc++.h>
#include <string.h>
#include <random>

using namespace std;

double Quadratic::f_quadratic(double x){

    return 0.3*x*x + 1.37*x + 5.43;
}

void Quadratic::create_quadratic_data(){

    size = 32;

    default_random_engine generator;
    random_device rd;
    mt19937 gen(rd());

    double xi = -5;
    double step = 10.0/double(size);
    int i = 0;

    while (xi < 5) {

        x.push_back(xi);
        normal_distribution<double> distribution( f_quadratic(x[i]) , 1.0 );
        y.push_back(distribution(generator));

        xi += step;
        i++;
    }
}

void Quadratic::read_quadratic_data(string data_file){

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

void Quadratic::fill_quadratic_system(){

    A.resize(3, vector<double>(3, 0));
    b.resize(3, 0);

    for (int i = 0; i < size; i++) {

        A[0][1] += x[i];
        A[0][2] += x[i]*x[i];
        A[1][0] += x[i];
        A[1][1] += x[i]*x[i];
        A[1][2] += x[i]*x[i]*x[i];
        A[2][0] += x[i]*x[i];
        A[2][1] += x[i]*x[i]*x[i];
        A[2][2] += x[i]*x[i]*x[i]*x[i];

        b[0] += y[i];
        b[1] += x[i]*y[i];
        b[2] += x[i]*x[i]*y[i];
    }

    A[0][0] = size;
}

void Quadratic::solve_quadratic(){

    coeffs.resize(3, 0);
    vector<vector<double> > A_inv(3, vector<double>(3, 0));
    double det = 0;

    det = A[0][0]*A[1][1]*A[2][2];
    det += A[0][1]*A[1][2]*A[2][0];
    det += A[0][2]*A[1][0]*A[2][1];
    det -= A[2][0]*A[1][1]*A[0][2];
    det -= A[2][1]*A[1][2]*A[0][0];
    det -= A[2][2]*A[1][0]*A[0][1];

    A_inv[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1])/det;
    A_inv[0][1] = (A[2][1]*A[0][2] - A[2][2]*A[0][1])/det;
    A_inv[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1])/det;
    A_inv[1][0] = (A[2][0]*A[1][2] - A[0][1]*A[2][2])/det;
    A_inv[1][1] = (A[0][0]*A[2][2] - A[2][0]*A[0][2])/det;
    A_inv[1][2] = (A[1][0]*A[0][2] - A[0][0]*A[1][2])/det;
    A_inv[2][0] = (A[1][0]*A[2][1] - A[2][0]*A[1][1])/det;
    A_inv[2][1] = (A[2][0]*A[0][1] - A[0][0]*A[2][1])/det;
    A_inv[2][2] = (A[0][0]*A[1][1] - A[1][0]*A[0][1])/det;

    for (int i = 0; i < 3; i++) {

        double sum = 0;

        for (int j = 0; j < 3; j++) {

            sum += A_inv[i][j]*b[j];
        }

        coeffs[i] = sum;
    }
}
