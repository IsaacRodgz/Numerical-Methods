#include <iostream>
#include<bits/stdc++.h>
#include "Linear.hpp"

using namespace std;

// g++ -std=c++11 Linear.cpp test.cpp -o test && ./test

void test_Quadratic(){

    /*
        Linear model
    */

    int n = 100;

    Linear lm(n);

    lm.create_linear_data();

    /*
    for (int i = 0; i < n; i++) {

        cout << lm.x[i] << " " << lm.y[i] << endl;
    }
    */

    lm.fill_linear_system();

    lm.solve_linear();

    cout << "\nRegression line: " << lm.coeffs[0] << "x + " << lm.coeffs[1] << "\n" << endl;

}

int main(int argc, char const *argv[]) {

    test_Quadratic();

    return 0;
}
