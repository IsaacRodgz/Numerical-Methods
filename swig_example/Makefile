all:

	swig -c++ -python fibonacci.i
	g++ -fpic -c fibonacci.cpp fibonacci_wrap.cxx -I/home/isaac/anaconda3/include/python3.6m
	gcc -shared fibonacci.o fibonacci_wrap.o -o _fibonacci.so -lstdc++

clean:

	rm -f *.o *.py *.cxx
