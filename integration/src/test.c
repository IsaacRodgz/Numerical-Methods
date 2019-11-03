#include<stdio.h>
#include<stdlib.h>
#include "integration.h"
#define PI 3.14159266

int main(int argc, char const *argv[]) {

    if (argc == 1) {

        printf("\nOne missing parameter: integral_rule_number\n");

        return 1;
    }

    switch (atoi(argv[1])) {
        case 1:
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 20 points: %e\n", trapezoidal_rule(sin_x, 0, 4*PI, 20));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 100 points: %e\n", trapezoidal_rule(sin_x, 0, 4*PI, 100));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 1000 points: %e\n\n", trapezoidal_rule(sin_x, 0, 4*PI, 1000));
            break;

        case 2:
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 20 points: %e\n", simpson3_rule(sin_x, 0, 4*PI, 20));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 100 points: %e\n", simpson3_rule(sin_x, 0, 4*PI, 100));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 1000 points: %e\n\n", simpson3_rule(sin_x, 0, 4*PI, 1000));
            break;

        case 3:
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 20 points: %e\n", simpson8_rule(sin_x, 0, 4*PI, 20));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 100 points: %e\n", simpson8_rule(sin_x, 0, 4*PI, 100));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 1000 points: %e\n\n", simpson8_rule(sin_x, 0, 4*PI, 1000));
            break;

        case 4:
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 20 points: %e\n", boole_rule(sin_x, 0, 4*PI, 20));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 100 points: %e\n", boole_rule(sin_x, 0, 4*PI, 100));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 1000 points: %e\n\n", boole_rule(sin_x, 0, 4*PI, 1000));
            break;

        case 6:
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 20 points: %e\n", weddle_rule(sin_x, 0, 4*PI, 20));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 100 points: %e\n", weddle_rule(sin_x, 0, 4*PI, 100));
            printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 1000 points: %e\n\n", weddle_rule(sin_x, 0, 4*PI, 1000));
            break;
    }

    return 0;
}
