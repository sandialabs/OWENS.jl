//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: sparse1.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "sparse1.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

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
