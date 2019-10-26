#include<stdio.h>
#include "interp2.h"

int main(int argc, char const *argv[]) {
    /*
    int n = 4;
    double x[] = {45, 50, 55, 60};
    double f[] = {0.7071, 0.766, 0.8192, 0.866};

    int m = 3;
    double p[] = {52, 56, 58};

    double y[3] = {0.0, 0.0, 0.0};

    gregory_forward_interp_t(n, x, f, m, p, y);


    int n = 5;
    double x[] = {1891, 1901, 1911, 1921, 1931};
    double f[] = {46, 66, 81, 93, 101};

    int m = 1;
    double p[] = {1925};

    double y[1] = {0.0};

    gregory_backward_interp_t(n, x, f, m, p, y);
    */

    int n = 5;
    double x[] = {0, 4, 8, 12, 16};
    double f[] = {14, 24, 32, 35, 40};

    int m = 1;
    double p[] = {9};

    double y[1] = {0.0};

    //gauss_forward_interp_t(n, x, f, m, p, y);
    //gauss_backward_interp_t(n, x, f, m, p, y);
    stirling_interp_t(n, x, f, m, p, y);

    printf("\n");
    for (int i = 0; i < m; i++) {
        printf("%.17g ", y[i]);
    }
    printf("\n\n");

    return 0;
}
