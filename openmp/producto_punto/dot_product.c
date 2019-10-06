#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char const *argv[]){

    printf("\nProducto punto de vectores\n--------------------\n");

    int num_hilos = atoi(argv[1]);
    omp_set_num_threads(num_hilos);
	printf("\nSe ejecuta con %i hilo(s)\n", num_hilos);

    int size = 100000000;

    double* a = malloc(size*sizeof(double));
	double* b = malloc(size*sizeof(double));
    double sum = 0;
    double sum2 = 0;

	double t_ini = omp_get_wtime();

    #pragma omp parallel for reduction(+:sum)
    for( int i = 0 ; i < size ; i++ ){
        sum += a[i]*b[i];
    }

	double t_fin = omp_get_wtime();
    
    for( int i = 0 ; i < size ; i++ ){
        sum2 += a[i]*b[i];
    }

    // Verifica que ambas soluciones son correctas
    if( sum != sum2 )
        printf("\nError in parallel operation\n\n");
    else
        printf("\nDuración de la operación: %f segundos con vectores de tamaño %d\n\n", t_fin - t_ini, size);

    free(a);
    free(b);

	return 0;
}
