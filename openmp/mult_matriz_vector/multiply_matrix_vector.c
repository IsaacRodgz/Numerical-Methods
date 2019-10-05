#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char const *argv[]){

    printf("\nMultiplicación de matriz por vector\n--------------------\n");

    int num_hilos = atoi(argv[1]);
    omp_set_num_threads(num_hilos);
	printf("\nSe ejecuta con %i hilo(s)\n", num_hilos);

    int size = 22000;

    double** A = malloc( size * sizeof *A );

    A[0] = malloc( (size * size) * sizeof **A );
    for (int i = 1; i < size; i++)
        A[i] = A[i-1] + size;

    double* a = malloc(size * sizeof a);
    double* b = malloc(size * sizeof b);

	double t_ini = omp_get_wtime();

    #pragma omp parallel for
    for( int i = 0 ; i < size ; i++ ){

        double sum = 0;

        for (int j = 0; j < size; j++) {

            sum += A[i][j] * a[j];
        }

        b[i] = sum;
    }

	double t_fin = omp_get_wtime();

    printf("\nDuración de la operación: %f segundos con vectores de tamaño %d\n\n", t_fin - t_ini, size);

    free(a);
    free(b);
    free(A[0]);
    free(A);

	return 0;
}
