#include <stdlib.h>
#include <stdio.h>
#include "plot.h"

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

void plotData( double (*f)(double), double xmin, double xmax, int nIntervals, double xSolution){

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

    tempDataFileName = "Data";
    tempSolutionFileName = "solution";
    gnuplotPipe = popen("gnuplot -persist","w");

    if (gnuplotPipe) {

        fprintf(gnuplotPipe,"set title \"Jacobi\"\n");
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
