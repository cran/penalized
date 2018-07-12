#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP penalized_CoxFitCpp(SEXP, SEXP, SEXP);
extern SEXP penalized_Lasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP penalized_Ridge(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP penalized_StepLasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"penalized_CoxFitCpp", (DL_FUNC) &penalized_CoxFitCpp, 3},
  {"penalized_Lasso",     (DL_FUNC) &penalized_Lasso,     9},
  {"penalized_Ridge",     (DL_FUNC) &penalized_Ridge,     9},
  {"penalized_StepLasso", (DL_FUNC) &penalized_StepLasso, 9},
  {NULL, NULL, 0}
};

void R_init_penalized(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
