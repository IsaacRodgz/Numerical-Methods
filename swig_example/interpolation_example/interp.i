%module interp1 %{

    #define SWIG_FILE_WITH_INIT
    #include "src/interp.h"

%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (int DIM1, double* IN_ARRAY1) {(int n1, double* xis), (int n2, double* fis), (int n3, double* pts)};
%apply (int DIM1, double* ARGOUT_ARRAY1) {(int n4, double* yis)};

%inline %{
    void newton_interp(int n1, double* xis, int n2, double* fis, int n3, double* pts, int n4, double* yis) {

        newton_interp_t(n1, xis, fis, n3, pts, yis);
    }
%}

/*%include "src/interp.h"*/
