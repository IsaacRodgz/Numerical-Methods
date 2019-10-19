#include<stdio.h>
#include "interp.h"

// http://www.personal.psu.edu/jjb23/web/htmls/sl455SP12/ch3/CH03_4B.pdf

int main(int argc, char const *argv[]) {
    /*
    int n = 3;
    double x[] = {1.3, 1.6, 1.9};
    double f[] = {0.620086, 0.4554022, 0.2818186};
    double fp[] = {-0.5220232, -0.5698959, -0.5811571};

    int m = 1;
    double p[] = {1.5};

    double y[1] = {0.0};
    */

    int n = 2;
    double x[] = {0.0, 1.570796};
    double f[] = {0.0, 1.0};
    double fp[] = {1.0, 0.0};

    int m = 10;
    double p[] = {-1.0, -0.77777778, -0.55555556, -0.33333333, -0.11111111, 0.11111111,  0.33333333,  0.55555556,  0.77777778,  1.0};

    double y[m];

    hermite_interp_t(n, x, f, fp, m, p, y);

    printf("\n");
    for (int i = 0; i < m; i++) {
        printf("%lf ", y[i]);
    }
    printf("\n\n");

    return 0;
}
