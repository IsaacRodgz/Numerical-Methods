#include<stdio.h>
#include<stdlib.h>
#include "../integration.h"
#define PI 3.14159266

int main(int argc, char const *argv[]) {

    //printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 20 points: %e\n", trapezoidal_rule(sin_x, 0, 4*PI, 20));
    //printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 100 points: %e\n", trapezoidal_rule(sin_x, 0, 4*PI, 100));
    //printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 1000 points: %e\n\n", trapezoidal_rule(sin_x, 0, 4*PI, 1000));

    printf("\nIntegral of 1/(1+x^2) in interval [-1, 1] with n = 1: %e\n", trapezoidal_rule(f1, -1, 1, 1));
    printf("\nIntegral of sqrt(1+cos^2(x)) in interval [1, 2] with n = 1: %e\n\n", trapezoidal_rule(f2, 1, 2, 1));

    return 0;
}
