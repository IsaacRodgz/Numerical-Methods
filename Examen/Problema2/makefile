gauss:

	gcc -lm -O2 -c GaussPivot_solver_test.c
	gcc -lm -O2 -c solve_matrix_direct.c
	gcc -lm -O2 -c matrix_struct.c
	gcc -lm -o runTest GaussPivot_solver_test.o solve_matrix_direct.o matrix_struct.o

doolittle:

	gcc -O2 -c Doolittle_solver_test.c -lm
	gcc -O2 -c solve_matrix_direct.c -lm
	gcc -O2 -c matrix_struct.c -lm
	gcc -o runTest Doolittle_solver_test.o solve_matrix_direct.o matrix_struct.o

cholesky:

	gcc -lm -O2 -c CholeskyModified_solver_test.c
	gcc -lm -O2 -c solve_matrix_direct.c
	gcc -lm -O2 -c matrix_struct.c
	gcc -lm -o runTest CholeskyModified_solver_test.o solve_matrix_direct.o matrix_struct.o

clean:

	rm -f runTest *.o
