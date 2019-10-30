#include <iostream>
#include<bits/stdc++.h>
#include <fstream>
#include <string>
#include <omp.h>

using namespace std;

void read_matrix(const string matrix_file){

    // Read file with matrix

    ifstream file(matrix_file);
    string line;

    // Ignore headers

    getline(file, line);

    // Read dimension

    getline(file, line);

    int m;
    int n;
    int nnz;

    stringstream ss(line);

    ss >> m >> n >> nnz;

    // Matrix elements

    vector<int> Mi_coo;
    vector<int> Mj_coo;
    vector<double> Mv_coo;

    vector<double> M(nnz, 0);
    vector<int> IM(m+1, 0);
    vector<int> JM(nnz, 0);

    /*
    Read elements in COO format
    */

    int i_coord, j_coord;
    double val;

    cout << m << " " << n << " " << nnz << "\n" << endl;

    for (int i = 0; i < nnz; i++) {

        file >> i_coord >> j_coord >> val;

        Mi_coo.push_back(i_coord-1);
        Mj_coo.push_back(j_coord-1);
        Mv_coo.push_back(val);
    }

    file.close();

    /*
    for (int i = 0; i < Mi_coo.size(); i++) {

        printf("M[%d][%d] = %lf\n", Mi_coo[i], Mj_coo[i], Mv_coo[i]);
    }
    cout << endl;
    */

    /*
    Construct CSR matrix from COO
    */

    //compute number of non-zero entries per row

    for (int i = 0; i < nnz; i++) {

        IM[Mi_coo[i]]++;
    }

    /*
    for (int i = 0; i <= m; i++) {

        printf("IM[%d] = %d\n", i, IM[i]);
    }
    cout << endl;
    */

    // Make cummulative sum to get index of row begin

    for (int i = 0, cum_sum = 0; i < m; i++) {

        int temp = IM[i];
        IM[i] = cum_sum;
        cum_sum += temp;
    }
    IM[m] = nnz;

    /*
    for (int i = 0; i <= m; i++) {

        printf("IM[%d] = %d\n", i, IM[i]);
    }
    cout << endl;
    */

    // Write Mj_coo and Mv_coo into JM and M respectively

    for (int i = 0; i < nnz; i++) {

        int row = Mi_coo[i];
        int dest = IM[row];

        //printf("row: %d, dest: %d\n", row, dest);

        JM[dest] = Mj_coo[i];
        M[dest] = Mv_coo[i];

        //printf("JM[%d] = %d, M[%d] = %lf\n", dest, Mj_coo[i], dest, Mv_coo[i]);

        IM[row]++;

        //printf("Update IM[%d] to %d\n\n", row, IM[row]);
    }

    for (int i = 0, last = 0; i <= m; i++) {

        int temp = IM[i];
        IM[i] = last;
        last = temp;
    }

    /*
    IM.push_back(0);

    for (int i = 0; i < m; i++) {

        for (int j = 0; j < n; j++) {

            if (A[i][j] != 0) {

                M.push_back(A[i][j]);
                JM.push_back(j);
                NNZ += 1;
            }
        }

        IM.push_back(NNZ);
    }
    */

    /*
    for (int i = 0; i <= m; i++) {

        printf("IM[%d] = %d\n", i, IM[i]);
    }

    cout << endl;

    for (int i = 0; i < nnz; i++) {

        printf("JM[%d] = %d\n", i, JM[i]);
    }

    cout << endl;

    for (int i = 0; i < nnz; i++) {

        printf("M[%d] = %lf\n", i, M[i]);
    }
    */

    int col, size, found;

    for (int i = 0; i < m; i++) {

        if (IM[i] < IM[i+1]){
            col = JM[IM[i]];
            size = IM[i+1] - IM[i];
            found = 0;
        }
        else
            continue;

        for (int j = col; found < size; j++) {
            if (col == j){
                printf("%lf ", M[IM[i]+found]);
                found++;
                col = JM[IM[i]+found];
            }
            else
                printf("%lf ", 0.0);
        }
        printf("\n");
    }

    printf("\n");
    int row = 2;
    for (int i = IM[row]; i < IM[row+1] ; i++) {

        printf("M[%d][%d] = %lf\n", row, JM[i], M[i]);
    }

    vector<double> x(n, 1.0);
    vector<double> y(m, 0.0);

    printf("\n");
    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {

        double y0 = 0.0;

        for (int j = IM[i]; j < IM[i+1]; j++) {

            y0 += M[j]*x[JM[j]];
        }

        y[i] = y0;
    }

    printf("\n");
    for (int i = 0; i < m; i++) {
        printf("%lf ", y[i]);
    }
    printf("\n");

}
