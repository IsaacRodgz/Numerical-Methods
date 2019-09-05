#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "plot.h"

double f1(double x){

    // f(x) = x^2

    return x*x;
}

double f2(double x){

    // f(x) = x^2 - 2

    return x*x - 2;
}

double f3(double x){

    // f(x) = sin(x)

    return sin(x);
}

double f4(double x){

    // f(x) = 1/x^2

    return 1/(x*x);
}

double f5(double x){

    // f(x) = x^3 + 3x^2 + 2x

    return x*x*x + 3*x*x + 2*x;
}

double f6(double x){

    // f(x) = 1/x^2

    return 1/(x-2);
}

void plotFunction( double (*f)(double), double xmin, double xmax, int nIntervals, double xSolution){

    double stepSize = (xmax - xmin) / nIntervals;

    double* xData = (double*) malloc((nIntervals)*sizeof(double));
    double* yData = (double*) malloc((nIntervals)*sizeof(double));

    for (int i = 0; i < nIntervals; i++) {

        xData[i] = xmin;
        yData[i] = f(xmin);
        xmin += stepSize;

    }

    FILE *gnuplotPipe, *tempDataFile, *tempSolutionFile;
    char *tempDataFileName;
    char *tempSolutionFileName;
    double x, y;

    tempDataFileName = "f(x)";
    tempSolutionFileName = "root";
    gnuplotPipe = popen("gnuplot -persist","w");

    if (gnuplotPipe) {

        fprintf(gnuplotPipe,"set title \"f(x)\"\n");
        fprintf(gnuplotPipe,"set xlabel \"x\"\n");
        fprintf(gnuplotPipe,"set ylabel \"f(x)\"\n");

        // Draw x axis
        fprintf(gnuplotPipe,"set arrow from %lf,0 to %lf,0 nohead\n", xData[0], xData[nIntervals-1]);

        fprintf(gnuplotPipe,"plot \"%s\" smooth csplines lc rgb \"blue\", \"%s\" with points lc rgb \"red\"\n", tempDataFileName, tempSolutionFileName);
        fflush(gnuplotPipe);

        // Write data: "x f(x)" to file tempSolutionFileName

        tempDataFile = fopen(tempDataFileName, "w");

        for (int i = 0; i <= nIntervals; i++) {

            x = xData[i];
            y = yData[i];

            fprintf(tempDataFile,"%lf %lf\n", x, y);
        }

        fclose(tempDataFile);

        // Write root: "x_solution 0" to file tempSolutionFileName

        tempSolutionFile = fopen(tempSolutionFileName, "w");
        fprintf(tempSolutionFile,"%lf %lf\n", xSolution, 0.0);
        fclose(tempSolutionFile);

        printf("press enter to continue...");
        getchar();

        remove(tempDataFileName);
        remove(tempSolutionFileName);
        fprintf(gnuplotPipe,"exit \n");

    } else {

        printf("gnuplot not found...");
    }
}

void plotData( double * data, int size ){

    FILE *gnuplotPipe, *tempDataFile, *tempSolutionFile;
    char *tempDataFileName;

    tempDataFileName = "Data";
    gnuplotPipe = popen("gnuplot -persist","w");

    if (gnuplotPipe) {

        fprintf(gnuplotPipe,"set title \"Data\"\n");
        fprintf(gnuplotPipe,"set xlabel \"x\"\n");
        fprintf(gnuplotPipe,"set yrange [0.5:-0.5]\n");

        // Draw x axis
        fprintf(gnuplotPipe,"set arrow from %lf,0 to %lf,0 nohead\n", data[0], data[size-1]);

        fprintf(gnuplotPipe,"plot \"%s\" with points lc rgb \"red\"\n", tempDataFileName);
        fflush(gnuplotPipe);

        // Write data: "x f(x)" to file tempSolutionFileName

        tempDataFile = fopen(tempDataFileName, "w");

        for (int i = 0; i < size; i++) {

            fprintf(tempDataFile,"%lf %lf\n", data[i], 0.0);
        }

        fclose(tempDataFile);

        printf("press enter to continue...");
        getchar();

        remove(tempDataFileName);
        fprintf(gnuplotPipe,"exit \n");

    } else {

        printf("gnuplot not found...");
    }
}
