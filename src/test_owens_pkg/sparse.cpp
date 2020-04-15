//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: sparse.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "sparse.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : emxArray_real_T *y_d
//                emxArray_int32_T *y_colidx
//                emxArray_int32_T *y_rowidx
// Return Type  : void
//
void b_sparse(emxArray_real_T *y_d, emxArray_int32_T *y_colidx, emxArray_int32_T
              *y_rowidx)
{
  int ctr;
  int col;
  int row;
  static const signed char x[144] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
  };

  ctr = y_d->size[0];
  y_d->size[0] = 12;
  emxEnsureCapacity_real_T(y_d, ctr);
  for (ctr = 0; ctr < 12; ctr++) {
    y_d->data[ctr] = 0.0;
  }

  ctr = y_colidx->size[0];
  y_colidx->size[0] = 13;
  emxEnsureCapacity_int32_T(y_colidx, ctr);
  for (ctr = 0; ctr < 13; ctr++) {
    y_colidx->data[ctr] = 0;
  }

  y_colidx->data[0] = 1;
  ctr = y_rowidx->size[0];
  y_rowidx->size[0] = 12;
  emxEnsureCapacity_int32_T(y_rowidx, ctr);
  for (ctr = 0; ctr < 12; ctr++) {
    y_rowidx->data[ctr] = 0;
  }

  y_rowidx->data[0] = 1;
  ctr = 0;
  for (col = 0; col < 12; col++) {
    for (row = 0; row < 12; row++) {
      if (x[row + 12 * col] != 0) {
        y_rowidx->data[ctr] = row + 1;
        y_d->data[ctr] = 1.0;
        ctr++;
      }
    }

    y_colidx->data[col + 1] = ctr + 1;
  }
}

//
// Arguments    : emxArray_real_T *y_d
//                emxArray_int32_T *y_colidx
//                emxArray_int32_T *y_rowidx
// Return Type  : void
//
void c_sparse(emxArray_real_T *y_d, emxArray_int32_T *y_colidx, emxArray_int32_T
              *y_rowidx)
{
  int ctr;
  int col;
  int row;
  static const signed char x[144] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
  };

  ctr = y_d->size[0];
  y_d->size[0] = 12;
  emxEnsureCapacity_real_T(y_d, ctr);
  for (ctr = 0; ctr < 12; ctr++) {
    y_d->data[ctr] = 0.0;
  }

  ctr = y_colidx->size[0];
  y_colidx->size[0] = 13;
  emxEnsureCapacity_int32_T(y_colidx, ctr);
  for (ctr = 0; ctr < 13; ctr++) {
    y_colidx->data[ctr] = 0;
  }

  y_colidx->data[0] = 1;
  ctr = y_rowidx->size[0];
  y_rowidx->size[0] = 12;
  emxEnsureCapacity_int32_T(y_rowidx, ctr);
  for (ctr = 0; ctr < 12; ctr++) {
    y_rowidx->data[ctr] = 0;
  }

  y_rowidx->data[0] = 1;
  ctr = 0;
  for (col = 0; col < 12; col++) {
    for (row = 0; row < 12; row++) {
      if (x[row + 12 * col] != 0) {
        y_rowidx->data[ctr] = row + 1;
        y_d->data[ctr] = 1.0;
        ctr++;
      }
    }

    y_colidx->data[col + 1] = ctr + 1;
  }
}

//
// Arguments    : const emxArray_real_T *varargin_1
//                emxArray_real_T *y_d
//                emxArray_int32_T *y_colidx
//                emxArray_int32_T *y_rowidx
//                int *y_m
//                int *y_n
// Return Type  : void
//
void d_sparse(const emxArray_real_T *varargin_1, emxArray_real_T *y_d,
              emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx, int *y_m,
              int *y_n)
{
  int mInt;
  int nInt;
  int numalloc;
  int row;
  int ctr;
  double xrc;
  mInt = varargin_1->size[0];
  nInt = varargin_1->size[1];
  numalloc = 0;
  row = varargin_1->size[0] * varargin_1->size[1];
  for (ctr = 0; ctr < row; ctr++) {
    if (varargin_1->data[ctr] != 0.0) {
      numalloc++;
    }
  }

  *y_m = varargin_1->size[0];
  *y_n = varargin_1->size[1];
  if (numalloc < 1) {
    numalloc = 1;
  }

  row = y_d->size[0];
  y_d->size[0] = numalloc;
  emxEnsureCapacity_real_T(y_d, row);
  for (row = 0; row < numalloc; row++) {
    y_d->data[row] = 0.0;
  }

  row = y_colidx->size[0];
  y_colidx->size[0] = varargin_1->size[1] + 1;
  emxEnsureCapacity_int32_T(y_colidx, row);
  ctr = varargin_1->size[1];
  for (row = 0; row <= ctr; row++) {
    y_colidx->data[row] = 0;
  }

  y_colidx->data[0] = 1;
  row = y_rowidx->size[0];
  y_rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(y_rowidx, row);
  for (row = 0; row < numalloc; row++) {
    y_rowidx->data[row] = 0;
  }

  y_rowidx->data[0] = 1;
  ctr = 0;
  for (numalloc = 0; numalloc < nInt; numalloc++) {
    for (row = 0; row < mInt; row++) {
      xrc = varargin_1->data[row + varargin_1->size[0] * numalloc];
      if (xrc != 0.0) {
        y_rowidx->data[ctr] = row + 1;
        y_d->data[ctr] = xrc;
        ctr++;
      }
    }

    y_colidx->data[numalloc + 1] = ctr + 1;
  }
}

//
// Arguments    : const double varargin_1[144]
//                emxArray_real_T *y_d
//                emxArray_int32_T *y_colidx
//                emxArray_int32_T *y_rowidx
// Return Type  : void
//
void sparse(const double varargin_1[144], emxArray_real_T *y_d, emxArray_int32_T
            *y_colidx, emxArray_int32_T *y_rowidx)
{
  int numalloc;
  int ctr;
  int row;
  double d;
  numalloc = 0;
  for (ctr = 0; ctr < 144; ctr++) {
    if (varargin_1[ctr] != 0.0) {
      numalloc++;
    }
  }

  if (numalloc < 1) {
    numalloc = 1;
  }

  ctr = y_d->size[0];
  y_d->size[0] = numalloc;
  emxEnsureCapacity_real_T(y_d, ctr);
  for (ctr = 0; ctr < numalloc; ctr++) {
    y_d->data[ctr] = 0.0;
  }

  ctr = y_colidx->size[0];
  y_colidx->size[0] = 13;
  emxEnsureCapacity_int32_T(y_colidx, ctr);
  for (ctr = 0; ctr < 13; ctr++) {
    y_colidx->data[ctr] = 0;
  }

  y_colidx->data[0] = 1;
  ctr = y_rowidx->size[0];
  y_rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(y_rowidx, ctr);
  for (ctr = 0; ctr < numalloc; ctr++) {
    y_rowidx->data[ctr] = 0;
  }

  y_rowidx->data[0] = 1;
  ctr = 0;
  for (numalloc = 0; numalloc < 12; numalloc++) {
    for (row = 0; row < 12; row++) {
      d = varargin_1[row + 12 * numalloc];
      if (d != 0.0) {
        y_rowidx->data[ctr] = row + 1;
        y_d->data[ctr] = d;
        ctr++;
      }
    }

    y_colidx->data[numalloc + 1] = ctr + 1;
  }
}

//
// File trailer for sparse.cpp
//
// [EOF]
//
