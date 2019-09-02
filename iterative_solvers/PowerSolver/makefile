all:

	gcc -c power_solver_test.c
	gcc -c ../solve_iterative.c
	gcc -c  ../matrix_struct.c
	gcc -c  ../plot.c
	gcc -o runTest power_solver_test.o solve_iterative.o matrix_struct.o plot.o -lm

clean:

	rm -f runTest power_solver_test.o solve_iterative.o matrix_struct.o plot.o
