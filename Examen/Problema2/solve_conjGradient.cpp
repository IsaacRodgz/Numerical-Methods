#include <iostream>
#include <math.h>
#include "conjGradient.h"

/*

Run command:

g++ -std=c++11 -O2 -fopenmp solve.cpp test.cpp -o runTest -O2 && ./runTest m.mtx

*/

using namespace std;

void test_solve(const string matrix_file, const string vector_file){

    vector<double> x = cGradientSolver( matrix_file, vector_file, 100000, 0.0000000001 );
	
    vector< vector<double> > A = read_matrix_full(matrix_file);

    vector<double> b = read_vector(vector_file);

    double error = 0;

    for(int i = 0; i < A.size(); i++){

        double sum = 0;

        for(int j = 0; j < A.size(); j++){

            sum += A[i][j]*x[j];
        }

        error += (sum - b[i])*(sum - b[i]);
    }

    printf("\nError: %lf\n", sqrt(error));
    
}

int main(int argc, char const *argv[]) {

    if ( argc == 1 ) {
        cout << "\nError: Faltan argumentos: matrix_file vector_file" << endl;
    }

    else if ( argc == 2 ) {
        cout << "\nError: Falta argumento: vector_file" << endl;
    }

    else{
        test_solve(argv[1], argv[2]);
    }

    return 0;
}
