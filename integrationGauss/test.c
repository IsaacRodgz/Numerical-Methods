#include<stdio.h>
#include<stdlib.h>
#include "integration.h"
#define PI 3.14159266

int main(int argc, char const *argv[]) {

    printf("\nIntegral of 1/(1+x^2) in interval [-1, 1] with 1 point: %e\n", one_point_gl(f1));
    printf("\nIntegral of 1/(1+x^2) in interval [-1, 1] with 2 points: %e\n", two_point_gl(f1));
    printf("\nIntegral of 1/(1+x^2) in interval [-1, 1] with 3 points: %e\n", three_point_gl(f1));
    printf("\nIntegral of 1/(1+x^2) in interval [-1, 1] with 4 points: %e\n\n", four_point_gl(f1));

    printf("\nIntegral of sqrt(1+cos^2(x)) in interval [1, 2] with 1 points: %e\n", one_point_gl(f2));
    printf("\nIntegral of sqrt(1+cos^2(x)) in interval [1, 2] with 2 points: %e\n", two_point_gl(f2));
    printf("\nIntegral of sqrt(1+cos^2(x)) in interval [1, 2] with 3 points: %e\n", three_point_gl(f2));
    printf("\nIntegral of sqrt(1+cos^2(x)) in interval [1, 2] with 4 points: %e\n\n", four_point_gl(f2));

    return 0;
}
