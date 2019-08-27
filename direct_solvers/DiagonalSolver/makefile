all:

	gcc -c diagonal_solver_test.c
	gcc -c ../solve_matrix_direct.c
	gcc -c  ../matrix_struct.c
	gcc -o runTest diagonal_solver_test.o solve_matrix_direct.o matrix_struct.o

clean:

	rm -f runTest diagonal_solver_test.o solve_matrix_direct.o matrix_struct.o
