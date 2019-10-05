#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char const *argv[]){

    printf("\nMultiplicación de vectores elemento a elemento\n--------------------\n");

    int num_hilos = atoi(argv[1]);
    omp_set_num_threads(num_hilos);
	printf("\nSe ejecuta con %i hilo(s)\n", num_hilos);

    int size = 100000000;

    double* a = malloc(size*sizeof(double));
	double* b = malloc(size*sizeof(double));
    double* c = malloc(size*sizeof(double));

	double t_ini = omp_get_wtime();

    #pragma omp parallel for
    for( int i = 0 ; i < size ; i++ ){
        c[i]=a[i]*b[i];
    }

	double t_fin = omp_get_wtime();

    printf("\nDuración de la operación: %f segundos con vectores de tamaño %d\n\n", t_fin - t_ini, size);

    free(a);
    free(b);
    free(c);

	return 0;
}
