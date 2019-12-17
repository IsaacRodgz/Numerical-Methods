#include <math.h>
#include <bits/stdc++.h>
#include <random>
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

    // L_sym laplacian

    for (int i = 0; i < graph.size(); i++) {

        for (int j = 0; j < graph.size(); j++) {

            if ( i == j ) {

                laplacian[i][j] = 1 - graph[i][j]/diag[i];
            }

            else {

                laplacian[i][j] = -graph[i][j]/(sqrt(diag[i]*diag[j]));
            }
        }
    }

}

// K-means algorithm

void Spectral::kMeansEig(int numIters, int numClusters){

    // Vector with the cluster assignment for each point

    clusterIDs.resize(eigenData.size(), 0);

    // Helper vect for not repeating centroids

    vector<int> chosen(eigenData.size(), 0);

    // Find random initial clusters

    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<std::mt19937::result_type> dist6(0,eigenData.size()-1);

    for (int i = 0; i < numClusters; i++) {

        while (true) {

            int cluster_index = dist6(rng);

            if (chosen[cluster_index] == 0) {

                chosen[cluster_index] = 1;

                vector<double> center(eigenData[0].size(), 0);

                for (int j = 0; j < center.size(); j++) {

                    center[j] = eigenData[cluster_index][j];
                }

                centers.push_back(center);

                break;
            }
        }
    }

    // Iterate and update centers

    for (int i = 0; i < numIters; i++) {

        //cout << "------------------------------------------------------------- "<< endl;
        //cout << "\n Iter " << i+1 << "\n" << endl;

        // Associate each point to the nearest center

        for (int j = 0; j < eigenData.size(); j++) {

            int newCenter = getNearestCenter(eigenData[j]);

            clusterIDs[j] = newCenter;
        }

        // Recompute centroid for each cluster

        vector<int> clusterSize(numClusters, 0);

        // Get number of points in each cluster

        for (int j = 0; j < eigenData.size(); j++) {

            int clustID = clusterIDs[j];

            clusterSize[clustID]++;
        }

        for (int j = 0; j < centers.size(); j++) {

            // Update only if cluster is not empty

            if ( clusterSize[j] > 0 ) {

                // Update each coordinate of the current centroid

                for (int k = 0; k < centers[j].size(); k++) {

                    // Iterate over all points in cluster

                    double new_val = 0;

                    for (int l = 0; l < eigenData.size(); l++) {

                        if ( clusterIDs[l] == j ) {

                            new_val += eigenData[l][k];
                        }
                    }

                    centers[j][k] = new_val*(1/static_cast<double>(clusterSize[j]));
                }
            }
        }

        //cout << "\n New clusters: \n" << endl;
        /*
        for (int k = 0; k < centers.size(); k++) {

            //cout << "\n Cluster " << k << " with " << clusterSize[k] << " points: \n" << endl;

            for (int l = 0; l < centers[k].size(); l++) {

                cout << centers[k][l] << ", ";
            }

            cout << endl;
        }
        */
    }
}

// Helper function for k-means algorithm

int Spectral::getNearestCenter(vector<double> p){

    double minDist = 1<<30;
    int currentCluster = 0;

    for (int i = 0; i < centers.size(); i++) {

        double dist = 0;

        for (int j = 0; j < p.size(); j++) {

            dist += pow(p[j]-centers[i][j], 2);
        }

        dist = sqrt(dist);

        if (dist < minDist) {

            minDist = dist;
            currentCluster = i;
        }
    }

    return currentCluster;
}

// Verify if A is diagonal

int Spectral::is_diagonal(vector< vector<double> > A, double epsilon){

    for(int i = 0; i < A.size(); i++){

        for (int j = 0; j < A[0].size(); j++) {

            if ( i != j ) {

                if( fabs( A[i][j] ) > epsilon )
                    return FALSE;
            }
        }
    }

    return TRUE;
}

void Spectral::computeEigen(int numEigen, int numIters, double epsilon){

    eigenVects.resize(laplacian.size(), vector<double>(laplacian.size(), 0));

    jacobiSolver(laplacian, eigenVects, numIters, epsilon);
}

void Spectral::transformData(int numEigen){

    vector<int> indexes(numEigen, 0);
    eigenVals_vect.resize(data.size(), 0);

    // Copy eigenvalues

    for (int i = 0; i < laplacian.size(); i++) {

        eigenVals_vect[i] = laplacian[i][i];
    }

    vector<double> eigenVals_copy = eigenVals_vect;

    sort(eigenVals_vect.begin(), eigenVals_vect.end());

    // Find indexes of eigenvectors to be used

    for (int i = 1; i < numEigen+1; i++) {

        int pos = -1;

        for (int j = 0; j < eigenVals_copy.size(); j++) {

            if ( eigenVals_vect[i] == eigenVals_copy[j] ) {

                indexes[i-1] = j;

                continue;
            }
        }
    }

    eigenData.resize(data.size(), vector<double>(numEigen, 0));

    for (int i = 0; i < data.size(); i++) {

        for (int j = 0; j < numEigen; j++) {

            eigenData[i][j] = eigenVects[i][indexes[j]];
        }
    }
    /*
    for (int i = 0; i < data.size(); i++) {

        for (int j = 0; j < data[0].size(); j++) {

            printf("%10lf ", data[i][j]);
        }

        if (eigenVects[i][indexes[0]] > 0) {

            printf("%10lf ", 1.0);
        }

        else {

            printf("%10lf ", 0.0);
        }


        cout << endl;
    }
    */
}

// Helper function to verify if A is symmetric

int Spectral::is_simetric(vector< vector<double> > &A){

    for(int i = 0; i < A.size(); i++){

        for(int j = 0; j < A[0].size(); j++){

            if(A[i][j] != A[j][i]){
                printf("%d, %d\n", i, j);
                printf("%e, %e\n", A[i][j], A[j][i]);
                return FALSE;
            }
        }
    }

    return TRUE;
}


double Spectral::vectNorm(vector<double> &vect){

    double sum = 0;

    for (int i = 0; i < vect.size(); i++) {
        sum += ( vect[i] * vect[i] );
    }

    if ( sum == 0.0 ) {
        fprintf(stderr, "\n[Error] Vector of norm 0 found\n\n");
    }

    return sqrt(sum);

}

void Spectral::maxOffDiagonal(vector< vector<double> > &A, int &p, int &q){

    double max = 0;

    for (int i = 0; i < A.size(); i++) {

        for (int j = 0; j < A[0].size(); j++) {

            if ( i != j) {

                if( fabs( A[i][j] ) > max){

                    max = fabs( A[i][j] );
                    p = i;
                    q = j;

                }
            }
        }
    }
}

void Spectral::jacobiSolver(vector< vector<double> > &A, vector< vector<double> > &F, int num_iters, double epsilon){

    if ( is_simetric(A) == FALSE ) {
        printf("Matrix is not symmetric, cannot apply algorithm\n");
        exit(-1);
    }

    for (int i = 0; i < F.size(); i++) {
        for (int j = 0; j < F[0].size(); j++) {
            if (i == j)
                F[i][j] = 1.0;
            else
                F[i][j] = 0.0;
        }
    }

    vector<double> FP(F.size(), 0);
    vector<double> FQ(F.size(), 0);

    int i;
    for (i = 0; i < num_iters; i++) {

        if ( (i+1) % 200 == 0 ) {

            cout << "\nIteration " << i+1 << endl;
        }

        // Get row, column position of max element of A off the diagonal
        int p;
        int q;

        maxOffDiagonal(A, p, q);

        if ( fabs( A[p][q] ) < epsilon ) {
            printf("\nSe llego a la convergencia en %d iteraciones\n", i+1);
            break;
        }

        // Copy columns p and q of matrix of eigenvectors F
        for (int k = 0; k < F.size(); k++) {
            FP[k] = F[k][p];
            FQ[k] = F[k][q];
        }

        // Calculate 2*a_ij / ( a_jj - a_ii )
        double w = (A[q][q] - A[p][p]) / (2*A[p][q]);
        double t = (1/(fabs(w) + sqrt(w*w + 1)))*( w > 0? 1 : -1 );

        double c = 1/sqrt(t*t + 1);
        double s = t*c;
        //double tau = s/(1+c);

        // Update matrix of eigenvectors F
        for (int k = 0; k < F.size(); k++) {
            F[k][p] = c*FP[k] - s*FQ[k];
            F[k][q] = s*FP[k] + c*FQ[k];
        }

        double temp_pq = A[p][q];
        A[p][q] = 0;
        A[q][p] = 0;

        double temp_pp = A[p][p];
        A[p][p] = (c*c)*temp_pp + (s*s)*A[q][q] - (2*c*s)*temp_pq;
        A[q][q] = (s*s)*temp_pp + (c*c)*A[q][q] + (2*c*s)*temp_pq;

        for (int j = 0; j < A.size(); j++) {
            double temp;
            if ( j != p & j != q ) {
                temp = A[j][p];
                A[j][p] = c*temp - s*A[j][q];
                A[j][q] = c*A[j][q] + s*temp;

                A[p][j] = A[j][p];
                A[q][j] = A[j][q];
            }
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

void Spectral::printEigen(){

    cout << "\n Eigenvects: \n" << endl;

    for (int i = 0; i < eigenVects.size(); i++) {

        for (int j = 0; j < eigenVects[0].size(); j++) {

            printf("%10lf ", eigenVects[i][j]);
        }

        cout << endl;
    }

    cout << "\n\n Eigenvals: \n" << endl;

    for (int i = 0; i < eigenVals.size(); i++) {

        for (int j = 0; j < eigenVals[0].size(); j++) {

            printf("%10lf ", eigenVals[i][j]);
        }

        cout << endl;
    }

    for (int i = 0; i < eigenVals_vect.size(); i++) {

        printf("%10lf ", eigenVals_vect[i]);
    }

    cout << endl;
}

void Spectral::printLaplacian(){

    cout << endl;

    for (int i = 0; i < laplacian.size(); i++) {

        for (int j = 0; j < laplacian[0].size(); j++) {

            printf("%10lf ", laplacian[i][j]);
        }

        cout << endl;
    }

    cout << endl;
}

void Spectral::printClusters(){

    ofstream file;
    file.open("output.txt");

    for (int i = 0; i < data.size(); i++) {

        for (int j = 0; j < data[0].size(); j++) {

            file << data[i][j] << " ";
        }

        file << " " << clusterIDs[i] << "\n";
    }

    file.close();
}
