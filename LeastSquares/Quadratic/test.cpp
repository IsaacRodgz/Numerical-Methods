#include <iostream>
#include<bits/stdc++.h>
#include "Quadratic.hpp"

using namespace std;

// g++ -std=c++11 Quadratic.cpp test.cpp -o test && ./test

void test_Quadratic(){

    Quadratic qm = Quadratic();

    qm.read_quadratic_data("quadratic_data.txt");

    qm.fill_quadratic_system();

    qm.solve_quadratic();

    cout << "\nRegression parabola: " << qm.coeffs[2] << "x^2 + " << qm.coeffs[1] << "x + " << qm.coeffs[0] << "\n" << endl;

}

int main(int argc, char const *argv[]) {

    test_Quadratic();

    return 0;
}
