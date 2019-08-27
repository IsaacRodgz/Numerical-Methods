all:

	gcc -c newton_solver_test.c
	gcc -c ../solve_iterative.c
	gcc -c ../matrix_struct.c
	gcc -c  ../plot.c
	gcc -o runTest newton_solver_test.o solve_iterative.o plot.o matrix_struct.o -lm

clean:

	rm -f runTest newton_solver_test.o solve_iterative.o plot.o matrix_struct.o
