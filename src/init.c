#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C++ source code.
*/

/* .C calls */
extern void radial(double *, double *, int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *);
extern void nonradial(double *, double *, int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *);

static const R_CMethodDef CEntries[] = {
    {"radial",    (DL_FUNC) &radial,    13},
    {"nonradial", (DL_FUNC) &nonradial, 13},
    {NULL, NULL, 0}
};

void R_init_npsf(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
