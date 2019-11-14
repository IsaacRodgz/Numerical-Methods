#include <iostream>
#include<bits/stdc++.h>
#include "LeastSquares.h"

using namespace std;

// g++ -std=c++11 LeastSquares.cpp test.cpp -o test && ./test

void test_LeastSquares(){

    int n = 100;

    LeastSquares lm(n);

    lm.create_data();

    /*
    for (int i = 0; i < n; i++) {

        cout << lm.x[i] << " " << lm.y[i] << endl;
    }
    */

    lm.fill_system();

    lm.solve();

    cout << "\nRegression line: " << lm.coeffs[0] << "x + " << lm.coeffs[1] << "\n" << endl;

}

int main(int argc, char const *argv[]) {

    test_LeastSquares();

    return 0;
}
