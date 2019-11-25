#include "FiniteElement.hpp"

#include <bits/stdc++.h>
#include <string.h>
#include <random>

using namespace std;

FiniteElement::FiniteElement(double lambda_p, int num_elements_p) : size(0), lambda(lambda_p),
    num_elements(num_elements_p) {}

double FiniteElement::f_linear(double x){

    return 1.37*x + 2.43;
}

void FiniteElement::create_linear_data(){

    size = 32;

    default_random_engine generator;
    random_device rd;
    mt19937 gen(rd());
    //uniform_real_distribution<> dis(-10, 10);

    double lenght = 10.0/size;

    for (double i = 0; i < 10; i+=lenght) {
        x.push_back(i);
        normal_distribution<double> distribution( f_linear(x[i]) , 0.3 );
        y.push_back(distribution(generator));

    }
}

void FiniteElement::read_linear_data(string data_file){

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

    size = x.size();
    file.close();
}

double FiniteElement::map(double t, double xmin, double xmax){

    return -1 + ((2)/(xmax-xmin))*(t-xmin);
}

double FiniteElement::n_i1(double r){

    return (1-r)/2;
}

double FiniteElement::n_i2(double r){

    return (1+r)/2;
}

void FiniteElement::fill_fe_system(){

    b.resize(num_elements, 0);
    A.resize(num_elements, vector<double>(num_elements, 0));

    // Calculate element indexes

    vector<int> indexes;
    int lenght = size/num_elements;

    indexes.push_back(0);

    for (int i = 1; i < num_elements; i++) {

        int sum = indexes[i-1] + lenght;
        indexes.push_back(sum);
    }

    // Calculate element system of equations

    // Elements from 0 to num_elements-1

    for (int i = 0; i < num_elements-1; i++) {

        for (int j = indexes[i]; j < indexes[i+1]; j++) {

            //cout << "TEST " << "i: " << i << "j: " << j << endl;

            double n1 = n_i1(map(x[j], x[indexes[i]], x[indexes[i+1]-1]));
            double n2 = n_i2(map(x[j], x[indexes[i]], x[indexes[i+1]-1]));

            b[i] += y[j]*n1;
            b[i+1] += y[j]*n2;

            A[i][i] += n1*n1;
            A[i][i+1] += n1*n2;
            A[i+1][i] += n2*n1;
            A[i+1][i+1] += n2*n2;

        }

        double len = x[indexes[i+1]]-x[indexes[i]];

        A[i][i] -= (lambda*len)/4;
        A[i][i+1] += (lambda*len)/4;
        A[i+1][i] += (lambda*len)/4;
        A[i+1][i+1] -= (lambda*len)/4;

    }

    for (int i = 0; i < indexes.size(); i++) {
        nodes.push_back(x[indexes[i]]);
    }
    nodes.push_back(x[indexes[indexes.size()-1]]);

}

void FiniteElement::solve_fem(){

        coeffs = b;

        // Temp array for p_k.T * A * p_k

        vector<double> temp( b.size(), 0 );

        // Calculate A * x

        for (int j = 0; j < A.size(); ++j){

            double sum = 0;

            for (int k = 0; k < A.size(); ++k){

                sum += A[j][k] * coeffs[k];

            }

            temp[j] = sum;
        }

        vector<double> r(b.size(), 0);

        // Initialize r_0 = b - A*x

        for (int i = 0; i < r.size(); i++) {
            r[i] = b[i] - temp[i];
        }

        vector<double> p(r.size(), 0);

        // Initialize p_0 = r_0

        for (int i = 0; i < p.size(); i++) {
            p[i] = r[i];
        }

        double alpha = 0;
        double beta = 0;

        double rNorm_old = 0;

        for (int j = 0; j < r.size(); ++j)
            rNorm_old += r[j] * r[j];

        for (int i = 0; i < num_elements; i++) {

            // Calculate A * p_k

            for (int j = 0; j < A.size(); ++j){

                double sum = 0;

                for (int k = 0; k < A.size(); ++k){

                    sum += A[j][k] * p[k];
                }

                temp[j] = sum;
            }

            // Calculate alpha_k = (p_k * r_k) / (p_k * A * p_k)

            double denom = 0;

            for (int j = 0; j < p.size(); ++j)
                denom += p[j] * temp[j];

            alpha = rNorm_old/denom;

            // Update x_k+1 = x_k + alpha*p_k

            for (int j = 0; j < coeffs.size(); j++) {
                coeffs[j] += alpha*p[j];
            }

            double rNorm_new = 0;

            // Update r_k+1 = r_k - alpha*A*p_k
            for (int j = 0; j < r.size(); j++) {
                r[j] -= alpha*temp[j];
                rNorm_new += r[j] * r[j];
            }

            if ( sqrt(rNorm_new) < 0.00000001 ) {
                //printf("\nAlgorithm converged after %d iterations\n\n", i+1);
                break;
            }

            beta = rNorm_new/rNorm_old;

            // Update p_k+1 = r_k+1 + beta*p_k

            for (int j = 0; j < p.size(); j++) {
                p[j] = r[j] + beta*p[j];
            }

            rNorm_old = rNorm_new;
        }

        vector<double> z(coeffs.size(), 0);

        for (int j = 0; j < A.size(); ++j){

            double sum = 0;

            for (int k = 0; k < A.size(); ++k){

                sum += A[j][k] * x[k];
            }

            z[j] = sum;
        }

        for (int i = 0; i < z.size(); i++) {
            z[i] = z[i] - b[i];
        }

        double sum = 0;

        for (int i = 0; i < z.size(); i++) {
            sum += z[i]*z[i];
        }
}
