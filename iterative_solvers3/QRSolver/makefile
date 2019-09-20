all:

	gcc -c qr_solver_test.c
	gcc -c ../solve_iterative2.c
	gcc -c ../solve_iterative3.c
	gcc -c  ../matrix_struct.c
	gcc -c ../solve_matrix_direct.c
	gcc -o runTest qr_solver_test.o solve_iterative3.o solve_iterative2.o solve_matrix_direct.o matrix_struct.o -lm -Wall

clean:

	rm -f runTest *.o
