#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char const *argv[]){

    printf("\nMultiplicación de matriz por matriz\n--------------------\n");

    int num_hilos = atoi(argv[1]);
    omp_set_num_threads(num_hilos);
	printf("\nSe ejecuta con %i hilo(s)\n", num_hilos);

    int size = 1000;

    double** A = malloc( size * sizeof *A );

    A[0] = malloc( (size * size) * sizeof **A );
    for (int i = 1; i < size; i++)
        A[i] = A[i-1] + size;

    double** B = malloc( size * sizeof *B );

    B[0] = malloc( (size * size) * sizeof **B );
    for (int i = 1; i < size; i++)
        B[i] = B[i-1] + size;

    double** C = malloc( size * sizeof *C );

    C[0] = malloc( (size * size) * sizeof **C );
    for (int i = 1; i < size; i++)
        C[i] = C[i-1] + size;

	double t_ini = omp_get_wtime();

    int k, l;

    #pragma omp parallel for private(k, l)
        for(int j = 0; j < size; j++){

            for(k = 0; k < size; k++){

                C[j][k] = 0;

                for(l = 0; l < size; l++){

                    C[j][k] += A[j][l] * B[l][k];
                }
            }
        }

	double t_fin = omp_get_wtime();

    printf("\nDuración de la operación: %f segundos con vectores de tamaño %d\n\n", t_fin - t_ini, size);

    free(A[0]);
    free(A);
    free(B[0]);
    free(B);
    free(C[0]);
    free(C);

	return 0;
}
