#include <stdio.h>
#include <stdlib.h>
#include "../interp.h"

void test_newton(const char *matrix_filename){

    double x[] = {0, 1, 2, 3, 5};
    double f[] = {-1, 1, 9, 29, 129};
    int n = 5;

    newton_interp(x, f, n);
}

void main(int argc, char* argv[]){
/*
    if(argc==1){

        printf("\nerror: Two missing arguments are required: /path/to/matrix, /path/to/vector\n\n");
    }

    if(argc==2){

        printf("\nerror: One missing argument is required: /path/to/vector\n");
    }

    if(argc==3){

        test_upper_triag(argv[1], argv[2]);
    }
*/

    test_newton("data.txt");

}
