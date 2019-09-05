#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../plot.h"

void test_solve_newton(char * funcName, char * xi){

    double xinit;

    sscanf(xi, "%lf", &xinit);

    double tolerance = 0.000000001;
    double epsilon = 0.000000001;
    int iterations = 1000;
    double x_solve;

    // Select function

    double (*func)(double);
    double (*funcP)(double);

    switch ( atoi(funcName) ) {
        case 1:
            func = &f1;
            funcP = &f1p;
            break;

        case 2:
            func = &f2;
            funcP = &f2p;
            break;

        case 3:
            func = &f3;
            funcP = &f3p;
            break;

        case 4:
            func = &f4;
            funcP = &f4p;
            break;

        case 5:
            func = &f5;
            funcP = &f5p;
            break;
    }

    x_solve = newtonSolver(func, funcP, xinit, tolerance, epsilon, iterations);

    printf("\n----------------------------------------------------\n");
    printf("\nRoot found by newton method, x_solve: %lf\n\n", x_solve);
    printf("----------------------------------------------------\n\n");

    printf("Plotting result...\n\n");

    double xRangeMin = x_solve - 2.0;
    double xRangeMax = x_solve + 2.0;
    int nIntervals = 1000;

    plotFunction(func, xRangeMin, xRangeMax, nIntervals, x_solve);

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Two missing arguments are required: [funcName] [initialValue]\n\n");
    }

    if(argc==2){

        printf("\nerror: One missing argument is required:[initialValue]\n\n");
    }

    else
        test_solve_newton(argv[1], argv[2]);
}
