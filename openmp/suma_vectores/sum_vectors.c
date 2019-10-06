#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int compare(double* a, double* b, int n){

    for(int i = 0; i < n; i++){

        if( a[i] != b[i] )
            return 0;
    }


    return 1;
}

int main(int argc, char const *argv[]){

    printf("\nSuma de vectores\n----------------\n");

    int num_hilos = atoi(argv[1]);
    omp_set_num_threads(num_hilos);
	printf("\nSe ejecuta con %i hilo(s)\n", num_hilos);

    int size = 100000000;

    double* a = malloc(size*sizeof(double));
	double* b = malloc(size*sizeof(double));
    double* c = malloc(size*sizeof(double));
    double* d = malloc(size*sizeof(double));

    for( int i = 0 ; i < size ; i++ ){
        d[i]=a[i]+b[i];
    }

	double t_ini = omp_get_wtime();

    #pragma omp parallel for
    for( int i = 0 ; i < size ; i++ ){
        c[i]=a[i]+b[i];
    }

	double t_fin = omp_get_wtime();

    // Verifica que ambas soluciones son correctas
    if( compare(c, d, size) == 0 )
        printf("\nError in parallel operation\n\n");
    else
        printf("\nDuración de la operación: %f segundos con vectores de tamaño %d\n\n", t_fin - t_ini, size);

    free(a);
    free(b);
    free(c);
    free(d);

	return 0;
}
