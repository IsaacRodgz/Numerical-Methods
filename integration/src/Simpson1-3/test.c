#include<stdio.h>
#include<stdlib.h>
#include "../integration.h"
#define PI 3.14159266

int main(int argc, char const *argv[]) {

    printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 20 points: %e\n", simpson3_rule(sin_x, 0, 4*PI, 20));
    printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 100 points: %e\n", simpson3_rule(sin_x, 0, 4*PI, 100));
    printf("\nIntegral of sin(x + 0.5) in interval [0, 4*pi] with 1000 points: %e\n\n", simpson3_rule(sin_x, 0, 4*PI, 1000));

    return 0;
}
