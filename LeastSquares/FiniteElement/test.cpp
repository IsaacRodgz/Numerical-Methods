#include <iostream>
#include<bits/stdc++.h>
#include "FiniteElement.hpp"

using namespace std;

// g++ -std=c++11 FiniteElement.cpp test.cpp -o test && ./test

void test_FE(){

    FiniteElement fem(0.01, 16);

    fem.read_linear_data("quadratic_data.txt");

    fem.fill_fe_system();

    fem.solve_fem();

    cout << "\nNodes and coefficients found:\n" << endl;

    for (int i = 0; i < fem.coeffs.size(); i++) {

        cout << fem.nodes[i] << " " << fem.coeffs[i] << endl;
    }
}

int main(int argc, char const *argv[]) {

    test_FE();

    return 0;
}
