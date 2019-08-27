all:

	gcc -c Jacobi_solver_test.c
	gcc -c ../solve_iterative.c
	gcc -c  ../matrix_struct.c
	gcc -o runTest Jacobi_solver_test.o solve_iterative.o matrix_struct.o -lm

clean:

	rm -f runTest Jacobi_solver_test.o solve_iterative.o matrix_struct.o
