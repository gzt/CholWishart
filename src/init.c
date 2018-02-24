#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP lmvgamma(SEXP, SEXP);
extern SEXP mvdigamma(SEXP, SEXP);
extern SEXP rCholWishart(SEXP, SEXP, SEXP);
extern SEXP rInvCholWishart(SEXP, SEXP, SEXP);
extern SEXP rInvWishart(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"lmvgamma",        (DL_FUNC) &lmvgamma,        2},
    {"mvdigamma",       (DL_FUNC) &mvdigamma,       2},
    {"rCholWishart",    (DL_FUNC) &rCholWishart,    3},
    {"rInvCholWishart", (DL_FUNC) &rInvCholWishart, 3},
    {"rInvWishart",     (DL_FUNC) &rInvWishart,     3},
    {NULL, NULL, 0}
};

void R_init_matrixdist(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
