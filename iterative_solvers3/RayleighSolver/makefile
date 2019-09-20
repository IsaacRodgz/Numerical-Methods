all:

	gcc -c rayleigh_solver_test.c
	gcc -c ../solve_iterative2.c
	gcc -c ../solve_iterative3.c
	gcc -c  ../matrix_struct.c
	gcc -c ../solve_matrix_direct.c
	gcc -o runTest rayleigh_solver_test.o solve_iterative3.o solve_iterative2.o matrix_struct.o solve_matrix_direct.o -lm -Wall

clean:

	rm -f runTest *.o
