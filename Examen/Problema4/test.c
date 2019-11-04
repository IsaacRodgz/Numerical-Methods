#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define PI 3.14159265358979323846
#define E 2.71828182845904523536

double func_true(double x){

    return -6.94226*pow(E, x) - 0.8*sin(x/2) - 0.4*cos(x/2);
}

double func(double x, double y){

    return y + sin(x/2);
}

double funcp(double x, double y){

    return y + sin(x/2) + 0.5*cos(x/2);
}

double hermite_integral(double y0, double y1, double y2, double yp0, double yp1, double h){

    return (h/15)*(7*y0+16*y1+7*y2) + ((h*h)/15)*(yp0-yp1);
}

int main(int argc, char const *argv[]) {

    double h = 0.001;
    int size = (2*PI)/h;

    double* y = malloc( size * sizeof *y );
    double* x = malloc( size * sizeof *x );
    double* z = malloc( size * sizeof *z );

    y[0] = 0.5;
    x[0] = -PI;

    y[1] = y[0] + h*func(x[0], y[0]);
    x[1] = x[0] + h;
    y[2] = y[1] + h*func(x[1], y[1]);
    x[2] = x[1] + h;

    for (int i = 3; i < size; i++) {

        y[i] = y[i-2] + hermite_integral( func(x[i-3], y[i-3]), func(x[i-2], y[i-2]), func(x[i-1], y[i-1]), funcp(x[i-3], y[i-3]), funcp(x[i-1] , y[i-1]), h);
        x[i] = x[i-1] + h;
    }

    for (int i = 0; i < size; i++) {
        printf("%lf %lf\n", x[i], y[i]);
    }

    double error = 0;

    for(int i = 0; i < size; i++){

        error += fabs(y[i]-func_true(x[i]));
    }

    //printf("Error of interpolating polynomial: %lf", error/size);

}
