solve:

	gcc -O2 -c Doolittle_solver_test.c
	gcc -O2 -c solve_matrix_direct.c
	gcc -O2 -c  matrix_struct.c
	gcc -O2 -o runTest Doolittle_solver_test.o solve_matrix_direct.o matrix_struct.o

inverse:

	gcc -O2 -c Inverse_solver_test.c
	gcc -O2 -c solve_matrix_direct.c
	gcc -O2 -c  matrix_struct.c
	gcc -O2 -o runTest Inverse_solver_test.o solve_matrix_direct.o matrix_struct.o

clean:

	rm -f runTest *.o
