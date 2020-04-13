//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: mtimes1.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//
#ifndef MTIMES1_H
#define MTIMES1_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void b_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
  *a_colidx, const emxArray_int32_T *a_rowidx, const double b[144], double c[144]);
extern void c_sparse_mtimes(const double a[144], const emxArray_real_T *b_d,
  const emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, double c
  [144]);
extern void d_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
  *a_colidx, const emxArray_int32_T *a_rowidx, const double b_data[], const int
  b_size[2], double c_data[], int c_size[2]);
extern void e_sparse_mtimes(const double a_data[], const emxArray_real_T *b_d,
  const emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, double c
  [144]);
extern void f_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
  *a_colidx, const emxArray_int32_T *a_rowidx, const double b_data[], double c
  [12]);
extern void g_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
  *a_colidx, const emxArray_int32_T *a_rowidx, int a_m, const emxArray_real_T
  *b_d, const emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, int
  b_n, emxArray_real_T *c_d, emxArray_int32_T *c_colidx, emxArray_int32_T
  *c_rowidx, int *c_m, int *c_n);
extern void sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
  *a_colidx, const emxArray_int32_T *a_rowidx, const emxArray_real_T *b_d, const
  emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx,
  coder_internal_sparse *c);

#endif

//
// File trailer for mtimes1.h
//
// [EOF]
//
