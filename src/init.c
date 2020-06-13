#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C++ source code.
*/

/* .C calls */
extern void radial(double *, double *, int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *);
extern void nonradial(double *, double *, int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *, double *, double *, int *, int *);
extern void gtre_ll(int *, double *, double *, int *, int *, double *, int *, double *, double *, double *, int *, double *, int *, double *);
extern void gtre_grad(int *, double *, double *, int *, int *, double *, int *, double *, double *, double *, int *, double *, int *, double *, double *);
extern void gtre(int *, double *, double *, int *, int *, double *, int *, double *, double *, double *, int *, double *, int *, int *, double *, double *, double *);
extern void HaltonSeq(int *, int *, int *, double *);
extern void Primes(int *, int *, double *);


static const R_CMethodDef CEntries[] = {
    {"radial",    (DL_FUNC) &radial,    13},
    {"nonradial", (DL_FUNC) &nonradial, 17},
    {"gtre_ll",   (DL_FUNC) &gtre_ll,   14},
    {"gtre_grad", (DL_FUNC) &gtre_grad, 15},
    {"gtre",      (DL_FUNC) &gtre,      17},
    {"HaltonSeq", (DL_FUNC) &HaltonSeq, 4},
    {"Primes",    (DL_FUNC) &Primes,    3},
    {NULL, NULL, 0}
};

void R_init_npsf(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
