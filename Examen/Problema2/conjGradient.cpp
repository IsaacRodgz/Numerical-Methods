#include <iostream>
#include<bits/stdc++.h>
#include <fstream>
#include <string>
#include <omp.h>
#define TRUE 1
#define FALSE 0

using namespace std;

int is_simetric(vector< vector<double> > A){

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

vector<double> read_vector(const string vector_file){

    // Read file with matrix

    ifstream file(vector_file);
    string line;

    // Read dimension

    getline(file, line);

    int m;
    int n;

    stringstream ss(line);

    ss >> m >> n;

    cout << "Read vector of size: " << m << "\n" << endl;

    double val;
    vector<double> b;

    for (int i = 0; i < m; i++) {

        file >> val;

        b.push_back(val);
    }

    file.close();

    return b;
}

vector< vector<double> > read_matrix(const string matrix_file){

    // Read file with matrix

    ifstream file(matrix_file);
    string line;

    // Read dimension

    getline(file, line);

    int m;
    int n;

    stringstream ss(line);

    ss >> m >> n;

    /*
    Read elements in COO format
    */

    cout << m << " " << n << "\n" << endl;

    double val;
    vector< vector<double> > A(m, vector<double>(n, 0));

    for (int i = 0; i < m*m; i++) {

        //cout << i/m << " " << i%m << endl;

        file >> val;

        A[i/m][i%m] = val;
    }

    file.close();

    return A;
}

vector< vector<double> > read_matrix_full(const string matrix_file){

    int nnz = 0;

    ifstream in(matrix_file);
    string unused;
    while ( getline(in, unused) )
        nnz++;
    in.close();
    nnz = nnz - 1;

    // Read file with matrix

    ifstream file(matrix_file);
    string line;

    // Read dimension

    getline(file, line);

    int m;
    int n;

    stringstream ss(line);

    ss >> m >> n;

    cout << m << " " << n << "\n" << endl;

    double val;
    int i, j;
    vector< vector<double> > A(m, vector<double>(n, 0));

    for (int k = 0; k < nnz; k++) {

        //cout << i/m << " " << i%m << endl;

        file >> i >> j >> val;

        //cout << i << " " << j << " " << val << endl;

        A[i-1][j-1] = val;
    }

    file.close();

    is_simetric(A);

    return A;
}

tuple< vector<double>, vector<int>, vector<int> > read_matrix_rala(const string matrix_file){

    int nnz = 0;

    ifstream in(matrix_file);
    string unused;
    while ( getline(in, unused) )
        nnz++;
    in.close();
    nnz = nnz - 1;

    // Read file with matrix

    ifstream file(matrix_file);
    string line;

    // Read dimension

    getline(file, line);

    int m;
    int n;

    stringstream ss(line);

    ss >> m >> n;

    double val;
    int i, j;
    vector< vector<double> > A(m, vector<double>(n, 0));

    for (int k = 0; k < nnz; k++) {

        //cout << i/m << " " << i%m << endl;

        file >> i >> j >> val;

        //cout << i << " " << j << " " << val << endl;

        A[i-1][j-1] = val;
    }

    file.close();

    if ( is_simetric(A) == FALSE ) {

        printf("\nA is not symmetric. System cannot be solved by conjugate gradient.\n\n");
        exit(-1);
    }

    // Matrix elements

    vector<int> Mi_coo;
    vector<int> Mj_coo;
    vector<double> Mv_coo;

    /*
    Read elements in COO format
    */

    for (int i = 0; i < A.size(); i++) {

        for (int j = 0; j < A[0].size(); j++) {

            if ( A[i][j] != 0 ) {

                Mi_coo.push_back(i);
                Mj_coo.push_back(j);
                Mv_coo.push_back(A[i][j]);
            }
        }
    }

    file.close();

    nnz = Mi_coo.size();

    cout << "\nRead matrix of size: (" << m << "," << n << ") with " << nnz << " non zero entries." << "\n" << endl;

    vector<double> M(nnz, 0);
    vector<int> IM(m+1, 0);
    vector<int> JM(nnz, 0);

    /*
    Construct CSR matrix from COO
    */

    //compute number of non-zero entries per row

    for (int i = 0; i < nnz; i++) {

        IM[Mi_coo[i]]++;
    }

    // Make cummulative sum to get index of row begin

    for (int i = 0, cum_sum = 0; i < m; i++) {

        int temp = IM[i];
        IM[i] = cum_sum;
        cum_sum += temp;
    }
    IM[m] = nnz;


    // Write Mj_coo and Mv_coo into JM and M respectively

    for (int i = 0; i < nnz; i++) {

        int row = Mi_coo[i];
        int dest = IM[row];

        JM[dest] = Mj_coo[i];
        M[dest] = Mv_coo[i];

        IM[row]++;

    }

    for (int i = 0, last = 0; i <= m; i++) {

        int temp = IM[i];
        IM[i] = last;
        last = temp;
    }

    return make_tuple(M, JM, IM);
}

vector<double> cGradientSolver(const string matrix_file, const string vector_file, int numIters, double epsilon){

    vector<double> A;
    vector<int> JA;
    vector<int> IA;

    tie(A, JA, IA) = read_matrix_rala(matrix_file);

    vector<double> b = read_vector(vector_file);

    vector<double> x = read_vector(vector_file);

    // Temp array for p_k.T * A * p_k

    vector<double> temp( b.size(), 0 );

    // Calculate A * x

    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int j = 0; j < b.size(); ++j){

        double sum = 0;

        for (int k = IA[j]; k < IA[j+1]; ++k){

            sum += A[k] * x[JA[k]];

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

    #pragma omp parallel for reduction(+:rNorm_old)
    for (int j = 0; j < r.size(); ++j)
        rNorm_old += r[j] * r[j];

    for (int i = 0; i < numIters; i++) {

        // Calculate A * p_k

        #pragma omp parallel for
        for (int j = 0; j < b.size(); ++j){

            double sum = 0;

            for (int k = IA[j]; k < IA[j+1]; ++k){

                sum += A[k] * p[JA[k]];
            }
            temp[j] = sum;
        }

        // Calculate alpha_k = (p_k * r_k) / (p_k * A * p_k)

        double denom = 0;
        #pragma omp parallel for reduction(+:denom)
        for (int j = 0; j < p.size(); ++j)
            denom += p[j] * temp[j];

        alpha = rNorm_old/denom;

        // Update x_k+1 = x_k + alpha*p_k

        for (int j = 0; j < x.size(); j++) {
            x[j] += alpha*p[j];
        }

        double rNorm_new = 0;

        // Update r_k+1 = r_k - alpha*A*p_k
        for (int j = 0; j < r.size(); j++) {
            r[j] -= alpha*temp[j];
            rNorm_new += r[j] * r[j];
        }

        if ( sqrt(rNorm_new) < epsilon ) {
            printf("\nAlgorithm converged after %d iterations\n\n", i+1);
            break;
        }

        beta = rNorm_new/rNorm_old;

        // Update p_k+1 = r_k+1 + beta*p_k

        for (int j = 0; j < p.size(); j++) {
            p[j] = r[j] + beta*p[j];
        }

        rNorm_old = rNorm_new;
    }

    vector<double> z(x.size(), 0);

    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int j = 0; j < z.size(); ++j){

        double sum = 0;

        for (int k = IA[j]; k < IA[j+1]; ++k){

            sum += A[k] * x[JA[k]];
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


    //printf("Error: %e\n", sqrt(sum));

    return x;
}


