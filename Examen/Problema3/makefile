make potencia:

	gcc -c -O2 k_power_solver_test.c
	gcc -c -O2 solve_iterative.c
	gcc -c -O2 matrix_struct.c
	gcc -o runTest k_power_solver_test.o solve_iterative.o matrix_struct.o -lm

make potencia_inversa:

	gcc -c -O2 k_inverse_power_solver_test.c
	gcc -c -O2 solve_iterative.c
	gcc -c -O2 matrix_struct.c
	gcc -o runTest k_inverse_power_solver_test.o solve_iterative.o matrix_struct.o -lm

clean:

	rm -f runTest *.o
