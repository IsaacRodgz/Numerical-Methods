#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../plot.h"

void test_solve_bisection(){

    double xmin = 1;
    double xmax = 2;
    double x_solve;

    x_solve = bisectionSolver(squareRoot, xmin, xmax, 0.0001);

    printf("\n----------------------------------------------------\n");
    printf("\nRoot found by bisection method, x_solve: %lf\n\n", x_solve);
    printf("----------------------------------------------------\n\n");

    printf("Plotting result...\n\n");

    double xRangeMin = -2.0;
    double xRangeMax = 2.0;
    int nIntervals = 1000;

    plotFunction(squareRoot, xRangeMin, xRangeMax, nIntervals, x_solve);

}

void main(int argc, char* argv[]){

        test_solve_bisection();
}
