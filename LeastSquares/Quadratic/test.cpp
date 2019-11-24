#include <iostream>
#include<bits/stdc++.h>
#include "Quadratic.hpp"

using namespace std;

// g++ -std=c++11 Quadratic.cpp test.cpp -o test && ./test

void test_Quadratic(){

    int n = 100;

    Quadratic qm(n);

    qm.create_quadratic_data();

    /*
    for (int i = 0; i < n; i++) {

        cout << qm.x[i] << " " << qm.y[i] << endl;
    }
    */

    qm.fill_quadratic_system();

    qm.solve_quadratic();

    cout << "\nRegression parabola: " << qm.coeffs[2] << "x^2 + " << qm.coeffs[1] << "x + " << qm.coeffs[0] << "\n" << endl;

}

int main(int argc, char const *argv[]) {

    test_Quadratic();

    return 0;
}
