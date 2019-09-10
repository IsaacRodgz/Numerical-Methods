all:

	gcc -c k_power_solver_test.c
	gcc -c ../solve_iterative.c
	gcc -c ../solve_matrix_direct.c
	gcc -c  ../matrix_struct.c
	gcc -o runTest k_power_solver_test.o solve_iterative.o solve_matrix_direct.o matrix_struct.o -lm

clean:

	rm -f runTest *.o
