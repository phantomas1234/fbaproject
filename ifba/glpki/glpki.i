%module glpki
%{
#include <glpk.h>
%}

/* Functions for creating C arrays in Python; see usage in 'example_refman.py' */
%include "carrays.i"
%array_class(int, intArray);
/* a = intArray(SIZE) -> "a" can be passed to C functions as int *, int [], ... */
%array_class(double, doubleArray);
/* a = doubleArray(SIZE) -> "a" can be passed to C functions as double *, double [], ... */

%include glpki.h