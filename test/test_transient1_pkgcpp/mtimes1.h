//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: mtimes1.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//
#ifndef MTIMES1_H
#define MTIMES1_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_transient1_types.h"

// Function Declarations
extern void b_sparse_mtimes(const double a[144], const emxArray_real_T *b_d,
  const emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, double c
  [144]);
extern void c_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
  *a_colidx, const emxArray_int32_T *a_rowidx, const double b_data[], double c
  [12]);
extern void sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
  *a_colidx, const emxArray_int32_T *a_rowidx, const double b[144], double c[144]);

#endif

//
// File trailer for mtimes1.h
//
// [EOF]
//
