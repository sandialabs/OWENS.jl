//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: str2double1.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//
#ifndef STR2DOUBLE1_H
#define STR2DOUBLE1_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern boolean_T b_copydigits(char s1[3], int *idx, char s, int *k, boolean_T
  allowpoint);
extern boolean_T b_copyexponent(char s1[3], int *idx, char s, int *k);
extern void b_readNonFinite(char s, int *k, boolean_T *b_finite, double *fv);
extern void b_readfloat(char s1[3], int *idx, char s, int *k, boolean_T
  allowimag, boolean_T *isimag, boolean_T *b_finite, double *nfv, boolean_T
  *foundsign, boolean_T *success);
extern void c_readfloat(emxArray_char_T *s1, int *idx, const emxArray_char_T *s,
  int *k, int n, boolean_T allowimag, boolean_T *isimag, boolean_T *b_finite,
  double *nfv, boolean_T *foundsign, boolean_T *success);
extern boolean_T copydigits(char s1[4], int *idx, const char s[2], int *k,
  boolean_T allowpoint);
extern boolean_T copyexponent(char s1[4], int *idx, const char s[2], int *k);
extern void d_readfloat(emxArray_char_T *s1, int *idx, const emxArray_char_T *s,
  int *k, int n, boolean_T allowimag, boolean_T *isimag, boolean_T *b_finite,
  double *nfv, boolean_T *foundsign, boolean_T *success);
extern void readNonFinite(const char s[2], int *k, boolean_T *b_finite, double
  *fv);
extern void readfloat(char s1[4], int *idx, const char s[2], int *k, boolean_T
                      allowimag, boolean_T *isimag, boolean_T *b_finite, double *
                      nfv, boolean_T *foundsign, boolean_T *success);

#endif

//
// File trailer for str2double1.h
//
// [EOF]
//
