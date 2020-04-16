//
// File: sparse1.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
//
#ifndef SPARSE1_H
#define SPARSE1_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void sparse_ctranspose(const emxArray_real_T *this_d, const
  emxArray_int32_T *this_colidx, const emxArray_int32_T *this_rowidx, int this_m,
  int this_n, emxArray_real_T *y_d, emxArray_int32_T *y_colidx, emxArray_int32_T
  *y_rowidx, int *y_m, int *y_n);
extern void sparse_spallocLike(int nzmax, emxArray_real_T *s_d, emxArray_int32_T
  *s_colidx, emxArray_int32_T *s_rowidx);
extern void sparse_sparse(int m, int n, int nzmaxval, emxArray_real_T *this_d,
  emxArray_int32_T *this_colidx, emxArray_int32_T *this_rowidx, int *this_m, int
  *this_n);

#endif

//
// File trailer for sparse1.h
//
// [EOF]
//
