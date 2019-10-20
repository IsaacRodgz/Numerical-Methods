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
    int n = 8;
    double x[] = {-2.0, -1.42857143, -0.85714286, -0.28571429, 0.28571429, 0.85714286, 1.42857143, 2.0};
    double f[] = {-0.90929743, -0.98990308, -0.75597537, -0.28184285, 0.28184285, 0.75597537, 0.98990308, 0.90929743};

    int m = 8;
    double p[] = {-1.0, -0.77777778, -0.55555556, -0.33333333, -0.11111111, 0.11111111,  0.33333333,  0.55555556};

    double y[m];

    newton_piecewise_interp_t(n, x, f, m, p, y);

    printf("\n");
    for (int i = 0; i < m; i++) {
        printf("%lf ", y[i]);
    }
    printf("\n\n");

    return 0;
}
