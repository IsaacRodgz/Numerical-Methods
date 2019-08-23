CC = gcc
CFLAGS = -Wall
OBJFILES = matrix_struct.o solve_matrix_direct.o test_methods.o
TARGET = test

all: $(TARGET)


$(TARGET) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES)

clean:
	rm -f $(OBJFILES) $(TARGET) *~

run:

	./test MATRICES/M_SIM.txt MATRICES/V_SIM.txt
