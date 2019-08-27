#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../plot.h"

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

    printf("Plotting result...\n\n");

    double xRangeMin = -2.0;
    double xRangeMax = 2.0;
    int nIntervals = 1000;

    plotData(squareRoot, xRangeMin, xRangeMax, nIntervals, x_solve);

}

void main(int argc, char* argv[]){

        test_solve_newton();
}
