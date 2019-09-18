all:

	gcc -c subspace_solver_test.c
	gcc -c ../solve_iterative3.c
	gcc -c ../solve_iterative2.c
	gcc -c  ../matrix_struct.c
	gcc -o runTest subspace_solver_test.o solve_iterative3.o solve_iterative2.o matrix_struct.o -lm -Wall

clean:

	rm -f runTest subspace_solver_test.o solve_iterative3.o solve_iterative2.o matrix_struct.o
