#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"

double polynomial4(double x){

    // f(x) = x^4 - 5x^2 - x - 3

    return (x*x*x*x) - 5*(x*x) - x - 3;
}

double polynomial4P(double x){

    // f(x) = 4x^3 - 10x - 1

    return 4*(x*x*x) - 10*(x) - 1;
}

double squareRoot(double x){

    // f(x) = x^2 - 2

    return x*x - 2;
}

double squareRootP(double x){

    // f'(x) = 2x

    return 2*x;
}

void test_solve_newton(){

    double xinit = 1;
    double tolerance = 0.000001;
    double epsilon = 0.0000000001;
    int iterations = 100;
    double x_solve;

    x_solve = newtonSolver(squareRoot, squareRootP, xinit, tolerance, epsilon, iterations);

    printf("\n----------------------------------------------------\n");
    printf("\nRoot found by newton method, x_solve: %lf\n\n", x_solve);
    printf("----------------------------------------------------\n\n");

}

void main(int argc, char* argv[]){

        test_solve_newton();
}
