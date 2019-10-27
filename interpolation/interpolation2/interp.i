%module interp2 %{

    #define SWIG_FILE_WITH_INIT
    #include "src/interp2.h"

%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (int DIM1, double* IN_ARRAY1) {(int n1, double* xis), (int n2, double* fis), (int n3, double* pts)};
%apply (int DIM1, double* ARGOUT_ARRAY1) {(int n4, double* yis)};

%inline %{
    void gregory_forward_interp(int n1, double* xis, int n2, double* fis, int n3, double* pts, int n4, double* yis) {

        gregory_forward_interp_t(n1, xis, fis, n3, pts, yis);
    }
%}

%inline %{
    void gregory_backward_interp(int n1, double* xis, int n2, double* fis, int n3, double* pts, int n4, double* yis) {

        gregory_backward_interp_t(n1, xis, fis, n3, pts, yis);
    }
%}

%inline %{
    void gauss_forward_interp(int n1, double* xis, int n2, double* fis, int n3, double* pts, int n4, double* yis) {

        gauss_forward_interp_t(n1, xis, fis, n3, pts, yis);
    }
%}

%inline %{
    void gauss_backward_interp(int n1, double* xis, int n2, double* fis, int n3, double* pts, int n4, double* yis) {

        gauss_backward_interp_t(n1, xis, fis, n3, pts, yis);
    }
%}

%inline %{
    void stirling_interp(int n1, double* xis, int n2, double* fis, int n3, double* pts, int n4, double* yis) {

        stirling_interp_t(n1, xis, fis, n3, pts, yis);
    }
%}
