all:

	gcc -c jacobi_solver_test.c
	gcc -c ../solve_iterative.c
	gcc -c ../solve_matrix_direct.c
	gcc -c  ../matrix_struct.c
	gcc -o runTest jacobi_solver_test.o solve_iterative.o solve_matrix_direct.o matrix_struct.o -lm -Wall

clean:

	rm -f runTest jacobi_solver_test.o solve_iterative.o solve_matrix_direct.o matrix_struct.o
