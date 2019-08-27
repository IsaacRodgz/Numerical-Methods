all:

	gcc -c gaussSeidel_solver_test.c -g
	gcc -c ../solve_iterative.c -g
	gcc -c  ../matrix_struct.c -g
	gcc -o runTest gaussSeidel_solver_test.o solve_iterative.o matrix_struct.o -lm -g

clean:

	rm -f runTest gaussSeidel_solver_test.o solve_iterative.o matrix_struct.o
