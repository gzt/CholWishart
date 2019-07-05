#ifndef R_CHOLWISHART_H
#define R_CHOLWISHART_H

#include <Rinternals.h> // For SEXP usage
#include <R_ext/Rdynload.h> // For R_GetCCallable

static R_INLINE SEXP  rCholWishart(SEXP ns, SEXP nuP, SEXP scal){
  // Initialize a static function pointer that persists between fn calls
  static SEXP (*fna)(SEXP, SEXP, SEXP) = NULL;
  // Only look it up once
  if (fna == NULL) {
    fna = (SEXP (*)(SEXP, SEXP, SEXP)) R_GetCCallable("CholWishart", "rCholWishart");
  }
  // Call the function
  return fna(ns, nuP, scal);
}

static R_INLINE SEXP  rInvCholWishart(SEXP ns, SEXP nuP, SEXP scal){
  // Initialize a static function pointer that persists between fn calls
  static SEXP (*fnb)(SEXP, SEXP, SEXP) = NULL;
  // Only look it up once
  if (fnb == NULL) {
    fnb = (SEXP (*)(SEXP, SEXP, SEXP)) R_GetCCallable("CholWishart", "rInvCholWishart");
  }
  // Call the function
  return fnb(ns, nuP, scal);
}

static R_INLINE SEXP  rInvWishart(SEXP ns, SEXP nuP, SEXP scal){
  // Initialize a static function pointer that persists between fn calls
  static SEXP (*fnc)(SEXP, SEXP, SEXP) = NULL;
  // Only look it up once
  if (fnc == NULL) {
    fnc = (SEXP (*)(SEXP, SEXP, SEXP)) R_GetCCallable("CholWishart", "rInvWishart");
  }
  // Call the function
  return fnc(ns, nuP, scal);
}

static R_INLINE SEXP lmvgamma(SEXP x, SEXP p){
  // Initialize a static function pointer that persists between fn calls
  static SEXP (*fnd)(SEXP, SEXP) = NULL;
  // Only look it up once
  if (fnd == NULL) {
    fnd = (SEXP (*)(SEXP, SEXP)) R_GetCCallable("CholWishart", "lmvgamma");
  }
  // Call the function
  return fnd(x, p);
}

static R_INLINE SEXP mvdigamma(SEXP x, SEXP p){
  // Initialize a static function pointer that persists between fn calls
  static SEXP (*fne)(SEXP, SEXP) = NULL;
  // Only look it up once
  if (fne == NULL) {
    fne = (SEXP (*)(SEXP, SEXP)) R_GetCCallable("CholWishart", "mvdigamma");
  }
  // Call the function
  return fne(x, p);
}

static R_INLINE SEXP  rPseudoWishart(SEXP ns, SEXP nuP, SEXP scal){
  // Initialize a static function pointer that persists between fn calls
  static SEXP (*fnf)(SEXP, SEXP, SEXP) = NULL;
  // Only look it up once
  if (fnf == NULL) {
    fnf = (SEXP (*)(SEXP, SEXP, SEXP)) R_GetCCallable("CholWishart", "rPseudoWishart");
  }
  // Call the function
  return fnf(ns, nuP, scal);
}

#endif // R_CHOLWISHART_H
