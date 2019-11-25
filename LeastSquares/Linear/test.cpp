#include <iostream>
#include<bits/stdc++.h>
#include "Linear.hpp"

using namespace std;

// g++ -std=c++11 Linear.cpp test.cpp -o test && ./test

void test_Lineal(){

    Linear lm = Linear();

    lm.read_linear_data("lineal_data.txt");

    lm.fill_linear_system();

    lm.solve_linear();

    cout << "\nRegression line: " << lm.coeffs[0] << "x + " << lm.coeffs[1] << "\n" << endl;

}

int main(int argc, char const *argv[]) {

    test_Lineal();

    return 0;
}
