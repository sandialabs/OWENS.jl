//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: sparse.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//
#ifndef SPARSE_H
#define SPARSE_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void b_sparse(emxArray_real_T *y_d, emxArray_int32_T *y_colidx,
                     emxArray_int32_T *y_rowidx);
extern void c_sparse(emxArray_real_T *y_d, emxArray_int32_T *y_colidx,
                     emxArray_int32_T *y_rowidx);
extern void d_sparse(const emxArray_real_T *varargin_1, emxArray_real_T *y_d,
                     emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx, int
                     *y_m, int *y_n);
extern void sparse(const double varargin_1[144], emxArray_real_T *y_d,
                   emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx);

#endif

//
// File trailer for sparse.h
//
// [EOF]
//
