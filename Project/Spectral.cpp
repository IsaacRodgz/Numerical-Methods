#include <math.h>
#include <bits/stdc++.h>
#define TRUE 1
#define FALSE 0

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

    laplacian.resize(graph.size(), vector<double>(graph.size()));

    for (int i = 0; i < graph.size(); i++) {

        for (int j = 0; j < graph.size(); j++) {

            if ( i == j ) {

                laplacian[i][j] = graph[i][j] - diag[i];
            }

            else {

                laplacian[i][j] = graph[i][j];
            }
        }
    }
}

// Verify if A is diagonal

int Spectral::is_diagonal(vector< vector<double> > A){

    for(int i = 0; i < A.size(); i++){

        for (int j = 0; j < A[0].size(); j++) {

            if ( i != j ) {

                if( fabs( A[i][j] ) > 0.00001 )
                    return FALSE;
            }
        }
    }

    return TRUE;
}

void Spectral::computeEigen(int numEigen, int numIters, double epsilon){

    eigenVects.resize(laplacian.size(), vector<double>(numEigen, 0));
    eigenVals.resize(numEigen, vector<double>(numEigen, 0));

    subspaceSolver(numIters, epsilon, numEigen);
}

void Spectral::subspaceSolver(int num_iters, double epsilon, int k){

    // eigenVects[n][p] - matrix of eigenvectors of A
    // eigenVals[p][p] - matrix of eigenvalues of A

    // Q[n][p]

    vector< vector<double> > Q(eigenVects.size(), vector<double>(eigenVects[0].size(), 0));

    // R[p][p]

    vector< vector<double> > R(eigenVals.size(), vector<double>(eigenVals[0].size(), 0));

    // Temp matrix for multiply( Q^T, A )

    vector< vector<double> > T(Q[0].size(), vector<double>(laplacian[0].size(), 0));

    // Initialize matrix of eigenvectors eigenVects

    for (size_t i = 0; i < eigenVects.size(); i++) {
        for (int j = 0; j < eigenVects[0].size(); j++) {
            eigenVects[i][j] = (double)rand()/RAND_MAX*2.0-1.0;
        }
    }

    for (int i = 0; i < num_iters; i++) {

        //cout << "Iter: " << i+1 << endl;

        // Iterate in subspace

        QRFactor(eigenVects, Q, R);

        // T = multiply( Q^T, A )

        int l, k;
        //#pragma omp parallel for private(l,k)
        for(int j = 0; j < T.size(); j++){
            for(k = 0; k < T[0].size(); k++){

                T[j][k] = 0;

                for(l = 0; l < laplacian.size(); l++){
                    T[j][k] += Q[l][j] * laplacian[l][k];
                }
            }
        }

        // eigenVals = multiply( T, Q );

        //#pragma omp parallel for private(l,k)
        for(int j = 0; j < eigenVals.size(); j++){
            for(k = 0; k < eigenVals[0].size(); k++){

                eigenVals[j][k] = 0;

                for(l = 0; l < Q.size(); l++){
                    eigenVals[j][k] += T[j][l] * Q[l][k];
                }
            }
        }

        // FI = multiply( A, Q );

        //#pragma omp parallel for private(l,k)
        for(int j = 0; j < eigenVects.size(); j++){
            for(k = 0; k < eigenVects[0].size(); k++){

                eigenVects[j][k] = 0;

                for(l = 0; l < Q.size(); l++){
                    eigenVects[j][k] += laplacian[j][l] * Q[l][k];
                }
            }
        }

        if (is_diagonal(eigenVals) == TRUE) {
            printf("\nSubspace method converged after %d iterations\n", i);
            return;
        }

    }

    printf("\nSubspace method did not converged after %d iterations. Returning last computations\n", num_iters);
}

void Spectral::QRFactor(vector< vector<double> > &A, vector< vector<double> > &Q, vector< vector<double> > &R){

    // Calculate r_00

    double norm = 0;

    for (int i = 0; i < A.size(); i++) {
        norm += A[i][0] * A[i][0];
    }

    norm = sqrt(norm);
    R[0][0] = norm;

    // Calculate q_0

    for (int i = 0; i < Q.size(); i++) {
        Q[i][0] = A[i][0] * (1/norm);
    }

    // Vector to form new q's
    vector<double> a_temp(A.size(), 0);

    for (int i = 1; i < A[0].size(); i++) {

        //cout << "QR iter: " << i+1 << endl;

        // Initialize a_temp to column A[i]
        for (int j = 0; j < A.size(); j++) {
            a_temp[j] = A[j][i];
        }

        for (int j = 0; j < i; j++) {

            double dot = 0;

            for (int k = 0; k < A.size(); k++) {
                //printf("Q[%d], A[%d]\n", Q->cols*k + j, A->cols*k + size);
                dot += Q[k][j] * A[k][i];
            }

            R[j][i] = dot;

            for (int k = 0; k < A.size(); k++) {
                a_temp[k] -= dot * Q[k][j];
            }
        }

        norm = 0;
        for (int j = 0; j < A.size(); j++) {
            norm += a_temp[j] * a_temp[j];
        }
        norm = sqrt(norm);

        // Calculate q_i

        for (int j = 0; j < Q.size(); j++) {
            Q[j][i] = a_temp[j] * (1/norm);
        }

        // Calculate r_ii

        R[i][i] = 0;

        for (int j = 0; j < A.size(); j++) {
            //printf("Q[%d], A[%d]\n", Q->cols*k + j, A->cols*k + size);
            R[i][i] += Q[j][i] * A[j][i];
        }
    }
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
