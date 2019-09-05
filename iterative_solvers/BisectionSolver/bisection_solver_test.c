#include <stdio.h>
#include <stdlib.h>
#include "../solve_iterative.h"
#include "../plot.h"

void test_solve_bisection(char * funcName, char * xi, char * xf){

    // Define interval to find solution

    double xmin;
    double xmax;

    sscanf(xi, "%lf", &xmin);
    sscanf(xf, "%lf", &xmax);

    // Select function

    double (*func)(double);

    switch ( atoi(funcName) ) {
        case 1:
            func = &f1;
            break;

        case 2:
            func = &f2;
            break;

        case 3:
            func = &f3;
            break;

        case 4:
            func = &f4;
            break;

        case 5:
            func = &f5;
            break;

        case 6:
            func = &f6;
            break;
    }

    double x_solve;
    int numIters = 5000;

    x_solve = bisectionSolver(func, xmin, xmax, 0.000001, numIters);

    printf("\n----------------------------------------------------\n");
    printf("\nRoot found by bisection method, x_solve: %lf\n\n", x_solve);
    printf("----------------------------------------------------\n\n");

    printf("Plotting result...\n\n");

    int nIntervals = 1000;

    plotFunction(func, xmin, xmax, nIntervals, x_solve);

}

void main(int argc, char* argv[]){

    if(argc==1){

        printf("\nerror: Three missing arguments are required: [funcName] [intervalLeft] [intervalRight]\n\n");
    }

    if(argc==2){

        printf("\nerror: Two missing arguments are required: [intervalLeft] [intervalRight]\n\n");
    }

    if(argc==2){

        printf("\nerror: One missing arguments are required: [intervalRight]\n\n");
    }
    else
        test_solve_bisection(argv[1], argv[2], argv[3]);
}
