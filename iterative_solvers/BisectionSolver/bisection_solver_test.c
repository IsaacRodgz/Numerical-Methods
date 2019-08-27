#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"

double polynomial4(double x){

    // f(x) = x^4 - 5x^2 - x - 3

    return (x*x*x*x) - 5*(x*x) - x - 3;
}

double squareRoot(double x){

    // f(x) = x^2 - 2

    return x*x - 2;
}

void test_solve_bisection(){

    double xmin = 1;
    double xmax = 2;
    double x_solve;

    x_solve = bisectionSolver(squareRoot, xmin, xmax, 0.0001);

    printf("\n----------------------------------------------------\n");
    printf("\nRoot found by bisection method, x_solve: %lf\n\n", x_solve);
    printf("----------------------------------------------------\n\n");

}

void main(int argc, char* argv[]){

        test_solve_bisection();
}
