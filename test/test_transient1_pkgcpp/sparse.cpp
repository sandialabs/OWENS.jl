//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: sparse.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "sparse.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "test_transient1_emxutil.h"
#include <string.h>

// Function Definitions

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
