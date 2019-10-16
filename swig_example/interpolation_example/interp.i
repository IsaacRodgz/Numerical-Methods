%module interp1 %{

    #define SWIG_FILE_WITH_INIT
    #include "src/interp.h"

%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (int DIM1, double* IN_ARRAY1) {(int n1, double* xis), (int n2, double* fis)};
%apply (int DIM1, double* ARGOUT_ARRAY1) {(int n3, double* coeffs)};

%inline %{
    void newton_interp_t(int n1, double* xis, int n2, double* fis, int n3, double* coeffs) {

        newton_interp(n1, xis, fis, coeffs);
    }
%}

/*%include "src/interp.h"*/
