/* file: hw.i */
%module swim
%{
/* Everything in this block will be copied in the wrapper file. We include the C header file necessary to compile the interface
*/
    #define SWIG_FILE_WITH_INIT
    #include "swim.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%apply (int DIM1, double* INPLACE_ARRAY1) {(int n1, double *a1), (int n2, double *a2), (int n3, double *a3)};

%include "swim.h"

/* list functions to be interfaced: */
