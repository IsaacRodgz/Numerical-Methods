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

// K-means algorithm

void Spectral::kMeans(int numIters, int numClusters){

    // Vector with the cluster assignment for each point

    clusterIDs.resize(data.size(), 0);

    // Helper vect for not repeating centroids

    vector<int> chosen(data.size(), 0);

    // Find random initial clusters

    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<std::mt19937::result_type> dist6(0,data.size()-1);

    for (int i = 0; i < numClusters; i++) {

        while (true) {

            int cluster_index = dist6(rng);

            if (chosen[cluster_index] == 0) {

                chosen[cluster_index] = 1;

                vector<double> center(data[0].size(), 0);

                for (int j = 0; j < center.size(); j++) {

                    center[j] = data[cluster_index][j];
                }

                centers.push_back(center);

                break;
            }
        }
    }

    // Iterate and update centers

    for (int i = 0; i < numIters; i++) {

        cout << "------------------------------------------------------------- "<< endl;
        cout << "\n Iter " << i+1 << "\n" << endl;

        // Associate each point to the nearest center

        for (int j = 0; j < data.size(); j++) {

            int newCenter = getNearestCenter(data[j]);

            clusterIDs[j] = newCenter;
        }

        // Recompute centroid for each cluster

        vector<int> clusterSize(numClusters, 0);

        // Get number of points in each cluster

        for (int j = 0; j < data.size(); j++) {

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

                    for (int l = 0; l < data.size(); l++) {

                        if ( clusterIDs[l] == j ) {

                            new_val += data[l][k];
                        }
                    }

                    centers[j][k] = new_val*(1/static_cast<double>(clusterSize[j]));
                }
            }
        }

        cout << "\n New clusters: \n" << endl;

        for (int k = 0; k < centers.size(); k++) {

            cout << "\n Cluster " << k << " with " << clusterSize[k] << " points: \n" << endl;

            for (int l = 0; l < centers[k].size(); l++) {

                cout << centers[k][l] << ", ";
            }

            cout << endl;
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

        cout << "------------------------------------------------------------- "<< endl;
        cout << "\n Iter " << i+1 << "\n" << endl;

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

        cout << "\n New clusters: \n" << endl;

        for (int k = 0; k < centers.size(); k++) {

            cout << "\n Cluster " << k << " with " << clusterSize[k] << " points: \n" << endl;

            for (int l = 0; l < centers[k].size(); l++) {

                cout << centers[k][l] << ", ";
            }

            cout << endl;
        }
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

// print cluster ID's

void Spectral::printClusters(){

    for (int i = 0; i < data.size(); i++) {

        for (int j = 0; j < data[0].size(); j++) {

            cout << data[i][j] << ", ";
        }

        cout << " ->  " << clusterIDs[i];

        cout << endl;
    }
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
    //eigenVals.resize(numEigen, vector<double>(numEigen, 0));
    //eigenVals_vect.resize(numEigen, 0);

    //subspaceSolver(numIters, epsilon, numEigen);

    //kInversePowerSolver(laplacian, eigenVects, eigenVals_vect, 10000, epsilon, numEigen);

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

    for (int i = 0; i < numEigen; i++) {

        int pos = -1;

        for (int j = 0; j < eigenVals_copy.size(); j++) {

            if ( eigenVals_vect[i] == eigenVals_copy[j] ) {

                indexes[i] = j;

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

        printf("%10lf ", eigenVects[i][indexes[0]]);

        cout << endl;
    }
    */
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

        if (is_diagonal(eigenVals, epsilon) == TRUE) {
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

// Solve linear system where A is a lower triangular matrix

vector<double> Spectral::solve_lower_triang(vector< vector<double> > &A, vector<double> &b, int fill_diag){

    if( A.size() != A[0].size() ){

        printf("\nNumber of rows and columns don't match. System cannot be solved\n\n");
        exit(-1);
    }

    if( is_diagonal(A, 0.000001) == FALSE ){

        printf("\nThere are diagonal elements equal to zero. System cannot be solved.\n\n");
        exit(-1);
    }

    double sum;
    int rows = A.size();
    int cols = A[0].size();

    // Vector of solution

    vector<double> x(rows, 0);

    // fill_diag = 1 changes diagional of A (A[i][i]) to 1's

    if(fill_diag == 1)
        x[0] = b[0];
    else
        x[0] = b[0]/A[0][0];

    // Calculate x[1] to x[n-1] with formula
    // x[i] = ( b[i] - Sum_(0 to i-1)_(A[i][j]*x[j]) ) / (A[i][i])

    for(int i = 1; i < rows; i++){

        sum = 0;

        for(int j = 0; j < i; j++){

            sum += A[i][j] * x[j];
        }

        if(fill_diag == 1)
            x[i] = b[i] - sum;
        else
            x[i] = (b[i] - sum) / ( A[i][i] ); // A[i][i] = A->data[ i*(rows+1) ]
    }

    return x;
}

vector<double> Spectral::solve_upper_triang(vector< vector<double> > &A, vector<double> &b, int fill_diag){

    if( A.size() != A[0].size() ){

        printf("\nNumber of rows and columns don't match. System cannot be solved\n\n");
        exit(-1);
    }

    if( is_diagonal(A, 0.000001) == FALSE ){

        printf("\nThere are diagonal elements equal to zero. System cannot be solved.\n\n");
        exit(-1);
    }

    double sum;
    int rows = A.size();
    int cols = A[0].size();

    // Vector of solution

    vector<double> x(rows, 0);

    // fill_diag = 1 changes diagional of A (A[i][i]) to 1's

    if(fill_diag == 1)
        x[rows-1] = b[rows-1];
    else
        x[rows-1] = b[rows-1] / A[rows-1][rows-1];

    // Calculate x[n-2] to x[0] with formula
    // x[i] = ( b[i] - Sum_(i+1 to n-1)_(A[i][j]*x[j]) ) / (A[i][i])

    for(int i = rows-2; i >= 0; i--){

        sum = 0;

        for(int j = i+1; j < rows; j++){

            sum += A[i][j] * x[j];
        }

        if(fill_diag == 1)
            x[i] = b[i] - sum;
        else
            x[i] = (b[i] - sum) / A[i][i]; // A[i][i] = A->data[ i*(n+1) ]
    }

    return x;
}

// Helper function to exchange rows from index: i_switch to index: index_set

void Spectral::set_pivot_row(vector< vector<double> > &A, vector<int> &pivot, int i_switch, int index_set){

    int cols = A[0].size();

    // Helper array to switch row and columns of A

    double buffer;

    // Switch rows in A

    for(int i = 0; i < cols; i++){

        buffer = A[i_switch][i];

        A[i_switch][i] = A[index_set][i];

        A[index_set][i] = buffer;
    }

    int temp = pivot[i_switch];
    pivot[i_switch] = pivot[index_set];
    pivot[index_set] = temp;

}

// Helper function to get the indexes of maximum element in current colum (in absolute value) of matrix A

int Spectral::max_pivot_index_column(vector< vector<double> > &A, int limit){

    int cols = A[0].size();
    int row = limit;
    double max = 0;

    for(int i = limit; i < cols; i++){

        if( fabs( A[i][limit] ) > max){

            max = fabs( A[i][limit] );
            row = i;

        }
    }

    return row;
}

// Helper function to factor matrix A into L*U with Doolittle and partial pivoting

void Spectral::factor_doolittle_pivoting(vector< vector<double> > &A, vector<int> &pivot){

    int cols = A[0].size();
    int current_pivot;

    for (int i = 0; i < cols; i++)
        pivot[i] = i;

    for (int i = 0; i < cols; i++) {

        current_pivot = max_pivot_index_column(A, i);

        if ( A[current_pivot][current_pivot] == 0.0 ) {
            printf("\nThere are diagonal elements equal to zero. System cannot be solved\n\n");
            exit(-1);
        }

        if ( current_pivot != i ) {

            set_pivot_row(A, pivot, current_pivot, i);
        }

        for (int j = i + 1; j < cols; j++) {

            A[j][i] = A[j][i] / A[i][i];

            for(int k = i + 1; k < cols; k++){

                A[j][k] -=  ( A[j][i] * A[i][k] );
            }
        }

    }
}

// Helper function to factor matrix A into L*U with Doolittle

void Spectral::factor_doolittle(vector< vector<double> > &A){

    int cols = A[0].size();

    // Alternate calculation of rows of U and columns of L

    for(int i = 0; i < cols; i++){

        // Calculate row i of U and save in A

        for(int j = i; j < cols; j++){

            // U[i][j] = A[i][j] - Sum_(k=0 to i-1){ L[i][k] * U[k][j] }

            for(int k = 0; k < i; k++){

                // A[i][j] -= A[i][k] * A[k][j]

                A[i][j] -= A[i][k] * A[k][j];
            }
        }

        // Calculate column i of L and save in A

        for(int j = i+1; j < cols; j++){

            // L[j][i] = A[j][i]/U[i][i] - Sum_(k=0 to i-1){ L[j][k] * U[k][i] }/U[i][i]

            // A[j][i] = A[j][i]/A[i][i]

            A[j][i] = A[j][i] / A[i][i];;

            for(int k = 0; k < i; k++){

                // A[j][i] -= (A[j][k] * A[k][i])/A[i][i]

                A[j][i] -= (( A[j][k] * A[k][i] ) / A[i][i]);
            }
        }
    }
}

vector<double> Spectral::solve_doolittle_pivoting(vector< vector<double> > &A, vector<double> &b, vector<int> &pivots){

    int rows = A.size();

    // Order b in x_solve

    vector<double> x_solve(rows, 0);

    for(int i = 0; i < rows; i++)
        x_solve[i] = b[pivots[i]];

    // Solve L*y = b where b is x_solve

    x_solve = solve_lower_triang(A, x_solve, 1);

    // Solve U*x = y where y is x_solve

    x_solve = solve_upper_triang(A, x_solve, 0);

    return x_solve;
}

vector<double> Spectral::solve_doolittle(vector< vector<double> > &A, vector<double> &b, int factor_flag){

    int rows = A.size();

    // Factor A as LU through Doolittle's method if factor_flag = 1

    if(factor_flag == 1)
        factor_doolittle(A);

    // Solve L*y = b

    vector<double> y(rows, 0);

    y = solve_lower_triang(A, b, 1);

    // Solve U*x = y

    vector<double> x_solve(rows, 0);

    x_solve = solve_upper_triang(A, y, 0);

    return x_solve;
}

// Find inverse of matrix A through Doolittle

vector< vector<double> > Spectral::solve_inverse(vector< vector<double> > &A){

    int cols = A[0].size();
    int rows = A.size();

    // Matrix to store inverse of A

    vector< vector<double> > inverse_A(rows, vector<double>(cols, 0));

    // Array of b. Takes columns of identity matrix

    vector<double> b(rows, 0);

    // Array to store solution of A*x = b, where x represents each column of A inverse

    vector<double> x(rows, 0);

    // Vector of pivots

    vector<int> pivots(rows, 0);

    // Factor A to LU through doolittle's algorithm

    factor_doolittle_pivoting(A, pivots);

    // Solve n linear sistems

    for(int i = 0; i < cols; i++){

        // Initialize b to i'th column of identity matrix

        for(int j = 0; j < cols; j++){

            if(j == i)
                b[j] = 1;
            else
                b[j] = 0;
        }

        // Solve for x with A already factored as LU

        x = solve_doolittle_pivoting(A, b, pivots);

        // Fill i'th column of A inverse with x

        for(int j = 0; j < cols; j++){

            inverse_A[j][i] = x[j];
        }
    }

    return inverse_A;
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

void Spectral::deflation(vector< vector<double> > &eigenVects, vector<double> &eigenVectInit, int currCol){

    vector<double> temp(eigenVectInit.size(), 0);

    for (int i = 0; i < currCol; i++) {

        double a = 0;

        for (int m = 0; m < eigenVects[0].size(); m++) {

            a += eigenVectInit[m] * eigenVects[i][m];
        }

        for (int m = 0; m < eigenVects[0].size(); m++) {

            temp[m] += a * eigenVects[i][m];
        }
    }

    for (int i = 0; i < eigenVectInit.size(); i++) {

        eigenVectInit[i] -= temp[i];
    }
}

void Spectral::kInversePowerSolver(vector< vector<double> > &A, vector< vector<double> > &eigenVects, vector<double> &eigenVals, int num_iters, double epsilon, int k){

    if ( is_simetric(A) == FALSE ) {
        printf("Matrix is not symmetric, cannot apply algorithm\n");
        exit(-1);
    }

    // Size of matrix and vectors
    int rows = A.size();

    // Variable to update dominant eigenvalue
    double lambdaNew;

    // Vector to calculate v_k
    vector<double> eigenVectNew;

    // Vector to calculate v_(k-1)
    vector<double> eigenVectOld(rows, 0);

    vector<int> pivots(rows, 0);
    factor_doolittle_pivoting(A, pivots);

    for (int s = 0; s < k; s++) {

        // Flag to indicate Convergence
        int converged = FALSE;

        // Initialize and Normalize initial vector v_0

        for (int i = 0; i < eigenVectOld.size(); i++) {
            eigenVectOld[i] = (double)rand()/RAND_MAX*2.0-1.0;
        }

        double norm;

        // Initialize lambda_0
        eigenVals[s] = 0;

        // Start iterations
        int i = 0;
        for (i = 0; i < num_iters; i++) {

            if (s>0)
                deflation(eigenVects, eigenVectOld, s);

            // Calculate v_(k) = w/norm(w)
            norm = vectNorm(eigenVectOld);
            for (int k = 0; k < rows; k++)
                eigenVectOld[k] = eigenVectOld[k] * (1/norm);

            // Calculate w =  A * v_(k-1)
            eigenVectNew = solve_doolittle_pivoting(A, eigenVectOld, pivots);

            // Calculate dominant eigenvalue and store in lambdaNew

            double accum2 = 0;
            double accum = 0;

            // v_k * (A*v_k)
            for (int j = 0; j < rows; j++) {
                accum += eigenVectNew[j] * eigenVectOld[j];
                accum2 += eigenVectNew[j] * eigenVectNew[j];
            }

            lambdaNew = accum/accum2;

            // Check for convergence and stop or update the eigenvalue

            if ( fabs( eigenVals[s] - lambdaNew ) < epsilon ) {
                eigenVals[s] = lambdaNew;
                converged = TRUE;
                break;
            }

            eigenVals[s] = lambdaNew;
            swap(eigenVectNew, eigenVectOld);
        }

        if (converged == TRUE) {
            printf("----------------------------------------------\n");
            printf("\nConverged for eigenvector_%d after %d iterations\n\n", s+1, i+1);

            for (int i = 0; i < eigenVects[0].size(); i++) {
                eigenVects[s][i] = eigenVectOld[i];
            }
        }
        else
            printf("\nMethod did not converge in given iterations. Returning last solution.\n\n");
    }
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

void Spectral::swap(vector<double> &A, vector<double> &B){

    for (int i = 0; i < A.size(); i++) {

        double temp = A[i];
        A[i] = B[i];
        B[i] = temp;
    }
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

        // Get row, column position of max element of A off the diagonal
        int p;
        int q;

        maxOffDiagonal(A, p, q);

        if ( fabs( A[p][q] ) < epsilon ) {
            printf("\nConverged at iteration: %d\n\n", i+1);
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

    printf("\nSe llego a la convergencia en %d iteraciones\n", i+1);
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

/*
void Spectral::printSpectral(){

    for (int i = 0; i < data.size(); i++) {

        for (int j = 0; j < data[0].size(); j++) {

            printf("%10lf ", data[i][j]);
        }

        printf("%10lf ", data[i][j]);

        cout << endl;
    }
}
*/
