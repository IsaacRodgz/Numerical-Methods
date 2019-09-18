#ifndef _matrix_struct_h
#define _matrix_struct_h

struct Matrix {
    int rows; // number of rows
    int cols; // number of columns
    double *data; // pointer to an array of rows*cols slots
};
typedef struct Matrix Matrix;

Matrix * read_matrix(const char *filename, int verbose);

void print_matrix(Matrix *A);

Matrix * multiply(Matrix *x, Matrix *y);

int equals(Matrix *x, Matrix *x_solve, double epsilon);

void copy(Matrix *x, Matrix *y);

int is_simetric(Matrix *A);

int is_diagonal(Matrix *A);

Matrix * transpose(Matrix *A);

#endif
