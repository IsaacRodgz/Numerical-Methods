/*
Author: Isaac Rodríguez Bribiesca
Date: 2019-08-12
Description: Pogram that solves system A*x = b of n equations with n variables,
through different methods.

Input:

    double *A : Pointer to matrix of nxn dimensions
    double *b : Pointer to array of n dimensions
    int n : Size of matrix A and array b

Output:

    double *x : Pointer to array of solution x

Methods included:

    - Diagonal matrix (solve_diagonal)
    - Lower triangular matrix (solve_lower_triang)
    - Upper triangular (solve_upper_triang)

*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "solve_matrix_direct.h"
#include "matrix_struct.h"
#define TRUE 1
#define FALSE 0

int is_diagonal(Matrix *A){

    int cols = A->rows;

    for(int i = 0; i < cols; i++){

        if( fabs( A->data[ (cols + 1)*i ] ) == 0.0 )
            return FALSE;
    }

    return TRUE;
}

Matrix * solve_diagonal(Matrix *A, Matrix *b){

    if( A->rows != A->cols ){

        printf("\nNumber of rows and columns don't match. System cannot be solved\n\n");
        exit(-1);
    }

    if( is_diagonal(A) != TRUE ){

        printf("\nThere are diagonal elements equal to zero. System cannot be solved.\n\n");
        exit(-1);
    }

    int rows = A->rows;

    // Vector of solution
    Matrix *x = malloc( sizeof( x ) );

    // Initialize Matrix struct
    x->rows = rows;
    x->cols = 1;
    x->data = malloc( x->rows*x->cols*sizeof( x->data ) );

    // Solution formula: x[i] = b[i] / A[i][i]

    for(int i = 0; i < rows; i++){

        x->data[i] = b->data[i] / ( A->data[ (rows + 1)*i ] );
    }

    return x;
}

Matrix * solve_lower_triang(Matrix *A, Matrix *b, int fill_diag){

    if( A->rows != A->cols ){

        printf("\nNumber of rows and columns don't match. System cannot be solved\n\n");
        exit(-1);
    }

    if( is_diagonal(A) != TRUE ){

        printf("\nThere are diagonal elements equal to zero. System cannot be solved.\n\n");
        exit(-1);
    }

    double sum;
    int rows = A->rows;
    int cols = A->cols;

    // Vector of solution
    Matrix *x = malloc( sizeof( x ) );

    // Initialize Matrix struct
    x->rows = rows;
    x->cols = 1;
    x->data = malloc( x->rows*x->cols*sizeof( x->data ) );

    // fill_diag = 1 changes diagional of A (A[i][i]) to 1's

    if(fill_diag == 1)
        x->data[0] = b->data[0];
    else
        x->data[0] = b->data[0]/A->data[0];

    // Calculate x[1] to x[n-1] with formula
    // x[i] = ( b[i] - Sum_(0 to i-1)_(A[i][j]*x[j]) ) / (A[i][i])

    for(int i = 1; i < rows; i++){

        sum = 0;

        for(int j = 0; j < i; j++){

            sum += A->data[ cols*i + j ] * x->data[ j ];
        }

        if(fill_diag == 1)
            x->data[i] = b->data[i] - sum;
        else
            x->data[i] = (b->data[i] - sum) / ( A->data[ i*(rows+1) ] ); // A[i][i] = A->data[ i*(rows+1) ]
    }

    return x;
}

Matrix * solve_upper_triang(Matrix *A, Matrix *b, int fill_diag){

    if( A->rows != A->cols ){

        printf("\nNumber of rows and columns don't match. System cannot be solved\n\n");
        exit(-1);
    }

    if( is_diagonal(A) != TRUE ){

        printf("\nThere are diagonal elements equal to zero. System cannot be solved.\n\n");
        exit(-1);
    }

    double sum;
    int rows = A->rows;

    // Vector of solution
    Matrix *x = malloc( sizeof( x ) );

    // Initialize Matrix struct
    x->rows = rows;
    x->cols = 1;
    x->data = malloc( x->rows*x->cols*sizeof( x->data ) );

    // fill_diag = 1 changes diagional of A (A[i][i]) to 1's

    if(fill_diag == 1)
        x->data[rows-1] = b->data[rows-1];
    else
        x->data[rows-1] = b->data[rows-1] / A->data[ rows*rows-1 ];

    // Calculate x[n-2] to x[0] with formula
    // x[i] = ( b[i] - Sum_(i+1 to n-1)_(A[i][j]*x[j]) ) / (A[i][i])

    for(int i = rows-2; i >= 0; i--){

        sum = 0;

        for(int j = i+1; j < rows; j++){

            sum += A->data[ rows*i + j ] * x->data[j];
        }

        if(fill_diag == 1)
            x->data[i] = b->data[i] - sum;
        else
            x->data[i] = (b->data[i] - sum) / A->data[ i*(rows+1) ]; // A[i][i] = A->data[ i*(n+1) ]
    }

    return x;
}

double diagonal_determinant(Matrix *A){

    int rows = A->rows;
    double determinant = 1;

    for(int i = 0; i < rows; i++){

        determinant *= A->data[ (rows + 1)*i ];
    }

    return determinant;
}

void set_pivot(Matrix *A, Matrix *b, int *index_order, int i_switch, int j_switch, int index_set){

    int cols = A->cols;

    // Helper array to switch row and columns of A

    double *buffer = malloc( cols*sizeof(*buffer) );

    // Helper variable to switch elements in b

    double b_temp;

    // Helper variable to switch index order of x_i

    int x_temp;

    if(i_switch != index_set){

        // Change in determinant sign

        index_order[cols] *= -1;

        // Switch elements in b

        b_temp = b->data[i_switch];
        b->data[i_switch] = b->data[index_set];
        b->data[index_set] = b_temp;

        // Switch rows in A

        for(int i = 0; i < cols; i++)
            buffer[i] = A->data[ i_switch*cols + i ];

        for(int i = 0; i < cols; i++)
            A->data[ i_switch*cols + i ] = A->data[ index_set*cols + i ];

        for(int i = 0; i < cols; i++)
            A->data[ index_set*cols + i ] = buffer[i];
    }

    if(j_switch != index_set){

        // Register order change in x_i

        x_temp = index_order[j_switch];
        index_order[j_switch] = index_order[index_set];
        index_order[index_set] = x_temp;

        // Switch columns in A

        for(int i = 0; i < cols; i++)
            buffer[i] = A->data[ i*cols + j_switch ];

        for(int i = 0; i < cols; i++)
            A->data[ i*cols + j_switch ] = A->data[ i*cols + index_set ];

        for(int i = 0; i < cols; i++)
            A->data[ i*cols + index_set ] = buffer[i];
    }

}

int * max_pivot_index(Matrix *A, int limit, int *max_ij){

    int cols = A->cols;

    double max = 0;

    for(int i = limit; i < cols; i++){

        for(int j = limit; j < cols; j++){

            if( fabs( A->data[ i*cols + j] ) > max){

                max = fabs( A->data[ i*cols + j] );
                max_ij[0] = i;
                max_ij[1] = j;

            }
        }
    }

    return max_ij;
}

void solve_gauss_elim(Matrix *A, Matrix *b){

    if( A->rows != A->cols ){

        printf("\nNumber of rows and columns don't match. System cannot be solved\n\n");
        exit(-1);
    }

    int cols = A->cols;

    // Reduce A to upper triangular

    for(int i = 0; i < cols-1; i++){

        if( A->data[ i*(cols+1) ] == 0.0 ){
            printf("\nThere are diagonal elements equal to zero. System cannot be solved\n\n");
            exit(-1);
        }

        for(int j = i+1; j < cols; j++){

            // A[j][i] = A[j][i] / A[i][i]

            A->data[ j*cols + i ] /=  A->data[ i*(cols+1) ];

            for(int k = i+1; k < cols; k++){

                // A[j][k] = A[j][k] - A[j][i] * A[i][k]

                A->data[ j*cols + k ] -= A->data[ j*cols + i ] * A->data[ i*cols + k ];
            }
        }
    }

    // Apply the same transformations to array b, stored in A[j][i]

    for(int i = 0; i < cols-1; i++){

        for(int j = i+1; j < cols; j++){

            b->data[j] -= A->data[ j*cols + i ] * b->data[i];
        }
    }
}

void solve_gauss_elim_pivot(Matrix *A, Matrix *b, int *index_order){

    if( A->rows != A->cols ){

        printf("\nNumber of rows and columns don't match. System cannot be solved\n\n");
        exit(-1);
    }

    int cols = A->cols;

    int *pivot_index = malloc( 2*sizeof( *pivot_index ) );

    // Reduce A to upper triangular

    for(int i = 0; i < cols-1; i++){

        // Find index pair (i,j) of max element from A

        max_pivot_index(A, i, pivot_index);
        /*
        printf("\n--------------------------------------------------\n");
        printf("\nmax[%d]: %f\n", cols*pivot_index[0] + pivot_index[1], A->data[cols*pivot_index[0] + pivot_index[1]]);// print max elem
        printf("\nbefore pivoting:\n");
        print_matrix(A);
        */
        set_pivot(A, b, index_order, pivot_index[0], pivot_index[1], i);
        /*
        printf("\nafter pivoting:\n");
        print_matrix(A);
        */
        if( A->data[ i*(cols+1) ] == 0.0 ){
            printf("\nMatrix is singular, cannot solve system\n");
            exit(-1);
        }

        for(int j = i+1; j < cols; j++){

            // A[j][i] = A[j][i] / A[i][i]

            A->data[ j*cols + i ] = A->data[ j*cols + i ] / A->data[ i*(cols+1) ];

            for(int k = i+1; k < cols; k++){

                // A[j][k] = A[j][k] - A[j][i] * A[i][k]

                A->data[ j*cols + k ] -= A->data[ j*cols + i ] * A->data[ i*cols + k ];
            }
        }
    }

    // Apply the same transformations to array b, stored in A[j][i]

    for(int i = 0; i < cols-1; i++){

        for(int j = i+1; j < cols; j++){

            b->data[j] -= A->data[ j*cols + i ] * b->data[i];
        }
    }

    free(pivot_index);
}

Matrix * solve_no_pivot(Matrix *A, Matrix *b){

    // Reduce matrix A to upper triangular through Gaussian elimination without pivoting (with the respective changes to b)

    solve_gauss_elim(A, b);

    // Perform back-substitution

    return solve_upper_triang(A, b, 0);
}

Matrix * solve_pivot(Matrix *A, Matrix *b){

    int cols = A->cols;
    int rows = A->rows;

    // Array to record reordering of solution x_i

    int * index_order = malloc( (rows + 1)*sizeof(*index_order) );

    // Matrix to save solution x with order of x_i saved in array index_order

    Matrix *x_solve = malloc( sizeof( x_solve ) );

    x_solve->rows = rows;
    x_solve->cols = 1;
    x_solve->data = malloc( x_solve->rows*x_solve->cols*sizeof( x_solve->data ) );

    // Matrix to return x in the correct order (x_0, x_1, ..., x_n)

    Matrix *x_solve_ordered = malloc( sizeof( x_solve_ordered ) );

    x_solve_ordered->rows = rows;
    x_solve_ordered->cols = 1;
    x_solve_ordered->data = malloc( x_solve_ordered->rows*x_solve_ordered->cols*sizeof( x_solve_ordered->data ) );

    // Initialize x indexes with 0, 1, ..., n to register change in pivot of columns. index_order[n] registers the sign for determinant

    for(int i = 0; i < cols; i++)
        index_order[i] = i;
    index_order[cols] = 1;

    // Reduce matrix A to upper triangular through Gaussian elimination with complete pivoting (with the respective changes to b)

    solve_gauss_elim_pivot(A, b, index_order);

    // Perform back-substitution

    x_solve = solve_upper_triang(A, b, 0);

    // Recover order of x_i through index_order and save it in x_solve_ordered

    for(int i = 0; i < cols; i++)
        x_solve_ordered->data[ index_order[i] ] = x_solve->data[i];

    free(index_order);
    free(x_solve);

    return x_solve_ordered;
}

void factor_doolittle(Matrix *A){

    int cols = A->cols;

    // Alternate calculation of rows of U and columns of L

    for(int i = 0; i < cols; i++){

        // Calculate row i of U and save in A

        for(int j = i; j < cols; j++){

            // U[i][j] = A[i][j] - Sum_(k=0 to i-1){ L[i][k] * U[k][j] }

            for(int k = 0; k < i; k++){

                // A[i][j] -= A[i][k] * A[k][j]

                A->data[ i*cols + j] -= A->data[ i*cols + k ] * A->data[ k*cols + j ];
            }
        }

        // Calculate column i of L and save in A

        for(int j = i+1; j < cols; j++){

            // L[j][i] = A[j][i]/U[i][i] - Sum_(k=0 to i-1){ L[j][k] * U[k][i] }/U[i][i]

            // A[j][i] = A[j][i]/A[i][i]

            A->data[ j*cols + i ] = A->data[ j*cols + i ] / A->data[ i*(cols+1) ];

            for(int k = 0; k < i; k++){

                // A[j][i] -= (A[j][k] * A[k][i])/A[i][i]

                A->data[ j*cols + i ] -= ( A->data[ j*cols + k] * A->data[ k*cols + i ] ) / A->data[ i*(cols+1) ];
            }
        }
    }
}

Matrix * solve_doolittle(Matrix *A, Matrix *b, int factor_flag){

    int rows = A->rows;

    // Factor A as LU through Doolittle's method if factor_flag = 1

    if(factor_flag == 1)
        factor_doolittle(A);

    // Solve L*y = b

    Matrix *y = malloc( sizeof( y ) );

    y->rows = rows;
    y->cols = 1;
    y->data = malloc( y->rows*y->cols*sizeof( y->data ) );

    y = solve_lower_triang(A, b, 1);

    // Solve U*x = y

    Matrix *x_solve = malloc( sizeof( x_solve ) );

    x_solve->rows = rows;
    x_solve->cols = 1;
    x_solve->data = malloc( x_solve->rows*x_solve->cols*sizeof( x_solve->data ) );

    x_solve = solve_upper_triang(A, y, 0);

    return x_solve;
}

Matrix * solve_inverse(Matrix *A){

    int cols = A->cols;
    int rows = A->rows;

    // Matrix to store inverse of A

    Matrix *inverse_A = malloc( sizeof( inverse_A ) );

    inverse_A->rows = rows;
    inverse_A->cols = cols;
    inverse_A->data = malloc( inverse_A->rows*inverse_A->cols*sizeof( inverse_A->data ) );

    // Array of b. Takes columns of identity matrix

    Matrix *b = malloc( sizeof( b ) );

    b->rows = rows;
    b->cols = 1;
    b->data = malloc( b->rows*b->cols*sizeof( b->data ) );

    // Array to store solution of A*x = b, where x represents each column of A inverse

    Matrix *x = malloc( sizeof( x ) );

    x->rows = rows;
    x->cols = 1;
    x->data = malloc( x->rows*x->cols*sizeof( x->data ) );

    // Factor A to LU through doolittle's algorithm

    factor_doolittle(A);

    // Solve n linear sistems

    for(int i = 0; i < cols; i++){

        // Initialize b to i'th column of identity matrix

        for(int j = 0; j < cols; j++){

            if(j == i)
                b->data[j] = 1;
            else
                b->data[j] = 0;
        }

        // Solve for x with A already factored as LU

        x = solve_doolittle(A, b, 0);

        // Fill i'th column of A inverse with x

        for(int j = 0; j < cols; j++){

            inverse_A->data[ j*cols + i] = x->data[j];
        }
    }

    return inverse_A;
}

void factor_cholesky_modified(Matrix *A){

    int cols = A->cols;

    for(int i = 0; i < cols; i++){

        // D[i][i] = A[i][i] - Sum_(k=1 to i-1){L[i][k]^2*D[k][k]}

        for(int k = 0; k < i; k++){

            A->data[ i*(cols+1) ] -= (A->data[ i*cols + k ] * A->data[ i*cols + k ] * A->data[ k*(cols+1) ]);
        }

        for(int j = i + 1; j < cols; j++){

            // L[j][i] = ( A[j][i] - Sum_(k=1 to i-1){L[j][k]*L[i][k]*D[k][k]} ) / D[i][i]

            A->data[ j*cols + i ] /= A->data[ i*(cols+1) ];

            for(int k = 0; k < i; k++){

                A->data[ j*cols + i ] -= ( ( A->data[ j*cols + k ] * A->data[ i*cols + k ] * A->data[ k*(cols+1) ] ) / A->data[ i*(cols+1) ] );
            }

            // Copy L[j][i] into its transpose L[i][j]

            A->data[ i*cols + j ] = A->data[ j*cols + i ];
        }
    }
}

Matrix * solve_cholesky_modified(Matrix *A, Matrix *b){

    int rows = A->rows;

    // Factor A as LDL^T through modified Cholesky method

    factor_cholesky_modified(A);

    // Solve for z in L*z = b

    Matrix *z = malloc( sizeof( z ) );

    z->rows = rows;
    z->cols = 1;
    z->data = malloc( z->rows*z->cols*sizeof( z->data ) );

    z = solve_lower_triang(A, b, 1);

    // // Solve for y in D*y = z

    Matrix *y = malloc( sizeof( y ) );

    y->rows = rows;
    y->cols = 1;
    y->data = malloc( y->rows*y->cols*sizeof( y->data ) );

    y = solve_diagonal(A, z);

    // // Solve for x in L^T*x = y

    Matrix *x_solve = malloc( sizeof( x_solve ) );

    x_solve->rows = rows;
    x_solve->cols = 1;
    x_solve->data = malloc( x_solve->rows*x_solve->cols*sizeof( x_solve->data ) );

    x_solve = solve_upper_triang(A, y, 1);

    return x_solve;
}
