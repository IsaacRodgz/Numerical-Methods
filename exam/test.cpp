#include <iostream>
#include "solve.h"

/*

Run command:

g++ -std=c++11 -fopenmp solve.cpp test.cpp -o runTest -O2 && ./runTest m.mtx

*/

using namespace std;

void test_solve(const string matrix_file){

    read_matrix(matrix_file);

}

int main(int argc, char const *argv[]) {

    if ( argc == 1 ) {
        cout << "\nError: Falta argumento: matrix_file" << endl;
    }

    else{
        test_solve(argv[1]);
    }

    return 0;
}
