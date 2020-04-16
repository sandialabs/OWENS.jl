//
// File: sparse1.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//

// Include Files
#include "sparse1.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_real_T *this_d
//                const emxArray_int32_T *this_colidx
//                const emxArray_int32_T *this_rowidx
//                int this_m
//                int this_n
//                emxArray_real_T *y_d
//                emxArray_int32_T *y_colidx
//                emxArray_int32_T *y_rowidx
//                int *y_m
//                int *y_n
// Return Type  : void
//
void sparse_ctranspose(const emxArray_real_T *this_d, const emxArray_int32_T
  *this_colidx, const emxArray_int32_T *this_rowidx, int this_m, int this_n,
  emxArray_real_T *y_d, emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx,
  int *y_m, int *y_n)
{
  int loop_ub;
  int idx;
  emxArray_int32_T *counts;
  int outridx_tmp;
  int outridx;
  sparse_sparse(this_n, this_m, this_colidx->data[this_colidx->size[0] - 1] - 1,
                y_d, y_colidx, y_rowidx, y_m, y_n);
  if ((this_m != 0) && (this_n != 0)) {
    loop_ub = y_colidx->size[0];
    for (idx = 0; idx < loop_ub; idx++) {
      y_colidx->data[idx] = 0;
    }

    idx = this_colidx->data[this_colidx->size[0] - 1];
    for (loop_ub = 0; loop_ub <= idx - 2; loop_ub++) {
      y_colidx->data[this_rowidx->data[loop_ub]]++;
    }

    y_colidx->data[0] = 1;
    idx = this_m + 1;
    for (loop_ub = 2; loop_ub <= idx; loop_ub++) {
      y_colidx->data[loop_ub - 1] += y_colidx->data[loop_ub - 2];
    }

    emxInit_int32_T(&counts, 1);
    idx = counts->size[0];
    counts->size[0] = this_m;
    emxEnsureCapacity_int32_T(counts, idx);
    for (idx = 0; idx < this_m; idx++) {
      counts->data[idx] = 0;
    }

    for (loop_ub = 0; loop_ub < this_n; loop_ub++) {
      for (idx = this_colidx->data[loop_ub] - 1; idx + 1 < this_colidx->
           data[loop_ub + 1]; idx++) {
        outridx_tmp = counts->data[this_rowidx->data[idx] - 1];
        outridx = (outridx_tmp + y_colidx->data[this_rowidx->data[idx] - 1]) - 1;
        y_d->data[outridx] = this_d->data[idx];
        y_rowidx->data[outridx] = loop_ub + 1;
        counts->data[this_rowidx->data[idx] - 1] = outridx_tmp + 1;
      }
    }

    emxFree_int32_T(&counts);
  }
}

//
// Arguments    : int nzmax
//                emxArray_real_T *s_d
//                emxArray_int32_T *s_colidx
//                emxArray_int32_T *s_rowidx
// Return Type  : void
//
void sparse_spallocLike(int nzmax, emxArray_real_T *s_d, emxArray_int32_T
  *s_colidx, emxArray_int32_T *s_rowidx)
{
  int numalloc;
  int i;
  if (nzmax >= 1) {
    numalloc = nzmax;
  } else {
    numalloc = 1;
  }

  i = s_d->size[0];
  s_d->size[0] = numalloc;
  emxEnsureCapacity_real_T(s_d, i);
  for (i = 0; i < numalloc; i++) {
    s_d->data[i] = 0.0;
  }

  i = s_colidx->size[0];
  s_colidx->size[0] = 13;
  emxEnsureCapacity_int32_T(s_colidx, i);
  s_colidx->data[0] = 1;
  i = s_rowidx->size[0];
  s_rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(s_rowidx, i);
  for (i = 0; i < numalloc; i++) {
    s_rowidx->data[i] = 0;
  }

  for (numalloc = 0; numalloc < 12; numalloc++) {
    s_colidx->data[numalloc + 1] = 1;
    s_colidx->data[numalloc] = 1;
  }

  s_colidx->data[12] = 1;
}

//
// Arguments    : int m
//                int n
//                int nzmaxval
//                emxArray_real_T *this_d
//                emxArray_int32_T *this_colidx
//                emxArray_int32_T *this_rowidx
//                int *this_m
//                int *this_n
// Return Type  : void
//
void sparse_sparse(int m, int n, int nzmaxval, emxArray_real_T *this_d,
                   emxArray_int32_T *this_colidx, emxArray_int32_T *this_rowidx,
                   int *this_m, int *this_n)
{
  int numalloc;
  int i;
  if (nzmaxval >= 1) {
    numalloc = nzmaxval;
  } else {
    numalloc = 1;
  }

  i = this_d->size[0];
  this_d->size[0] = numalloc;
  emxEnsureCapacity_real_T(this_d, i);
  for (i = 0; i < numalloc; i++) {
    this_d->data[i] = 0.0;
  }

  i = this_colidx->size[0];
  this_colidx->size[0] = n + 1;
  emxEnsureCapacity_int32_T(this_colidx, i);
  this_colidx->data[0] = 1;
  i = this_rowidx->size[0];
  this_rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(this_rowidx, i);
  for (i = 0; i < numalloc; i++) {
    this_rowidx->data[i] = 0;
  }

  for (numalloc = 0; numalloc < n; numalloc++) {
    this_colidx->data[numalloc + 1] = 1;
  }

  i = this_colidx->size[0];
  for (numalloc = 0; numalloc <= i - 2; numalloc++) {
    this_colidx->data[numalloc] = 1;
  }

  this_colidx->data[this_colidx->size[0] - 1] = 1;
  *this_m = m;
  *this_n = n;
}

//
// File trailer for sparse1.cpp
//
// [EOF]
//
