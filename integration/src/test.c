#include<stdio.h>
#include "integration.h"
#define PI 3.14159266

int main(int argc, char const *argv[]) {

    //printf("\nIntegral sin(x): %lf\n", trapezoidal_rule(sin_x, 0, 4*PI, 20));
    //printf("\nIntegral sin^2(x): %lf\n\n", trapezoidal_rule(sin_sqr_x, 0, PI, 6));

    //printf("\nIntegral sin(x): %lf\n", simpson3_rule(sin_x, 0, 4*PI, 20));
    //printf("\nIntegral sin^2(x): %lf\n\n", simpson3_rule(sin_sqr_x, 0, PI, 6));

    //printf("\nIntegral sin(x): %lf\n", simpson8_rule(sin_x, 0, 4*PI, 20));
    //printf("\nIntegral sin^2(x): %lf\n\n", simpson8_rule(sin_sqr_x, 0, PI, 6));

    //printf("\nIntegral sin(x): %lf\n", boole_rule(sin_x, 0, 4*PI, 20));
    //printf("\nIntegral sin^2(x): %lf\n\n", boole_rule(sin_sqr_x, 0, PI, 6));

    printf("\nIntegral sin(x): %lf\n", weddle_rule(sin_x, 0, 4*PI, 20));
    printf("\nIntegral sin^2(x): %lf\n\n", weddle_rule(sin_sqr_x, 0, PI, 6));

    return 0;
}
