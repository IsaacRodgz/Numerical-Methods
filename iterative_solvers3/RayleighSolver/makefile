all:

	gcc -c rayleigh_solver_test.c
	gcc -c ../solve_iterative3.c
	gcc -c  ../matrix_struct.c
	gcc -o runTest rayleigh_solver_test.o solve_iterative3.o matrix_struct.o -lm -Wall

clean:

	rm -f runTest rayleigh_solver_test.o solve_iterative3.o matrix_struct.o
