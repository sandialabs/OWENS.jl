//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: sort.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "sort.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include "sortLE.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : emxArray_creal_T *x
//                emxArray_int32_T *idx
// Return Type  : void
//
void b_sort(emxArray_creal_T *x, emxArray_int32_T *idx)
{
  int dim;
  emxArray_creal_T *vwork;
  int i;
  int vlen;
  int i2;
  int vstride;
  int k;
  emxArray_int32_T *iidx;
  emxArray_int32_T *iwork;
  emxArray_creal_T *xwork;
  int j;
  int n;
  int b_j;
  int pEnd;
  int p;
  int q;
  int qEnd;
  int kEnd;
  dim = 0;
  if (x->size[0] != 1) {
    dim = -1;
  }

  emxInit_creal_T(&vwork, 1);
  if (dim + 2 <= 1) {
    i = x->size[0];
  } else {
    i = 1;
  }

  vlen = i - 1;
  i2 = vwork->size[0];
  vwork->size[0] = i;
  emxEnsureCapacity_creal_T(vwork, i2);
  i = idx->size[0];
  idx->size[0] = x->size[0];
  emxEnsureCapacity_int32_T(idx, i);
  vstride = 1;
  for (k = 0; k <= dim; k++) {
    vstride *= x->size[0];
  }

  emxInit_int32_T(&iidx, 1);
  emxInit_int32_T(&iwork, 1);
  emxInit_creal_T(&xwork, 1);
  for (j = 0; j < vstride; j++) {
    for (k = 0; k <= vlen; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    i = vwork->size[0];
    n = vwork->size[0] + 1;
    i2 = iidx->size[0];
    iidx->size[0] = vwork->size[0];
    emxEnsureCapacity_int32_T(iidx, i2);
    dim = vwork->size[0];
    for (i2 = 0; i2 < dim; i2++) {
      iidx->data[i2] = 0;
    }

    if (vwork->size[0] != 0) {
      i2 = iwork->size[0];
      iwork->size[0] = vwork->size[0];
      emxEnsureCapacity_int32_T(iwork, i2);
      i2 = vwork->size[0] - 1;
      for (k = 1; k <= i2; k += 2) {
        if (sortLE(vwork, k, k + 1)) {
          iidx->data[k - 1] = k;
          iidx->data[k] = k + 1;
        } else {
          iidx->data[k - 1] = k + 1;
          iidx->data[k] = k;
        }
      }

      if ((vwork->size[0] & 1) != 0) {
        iidx->data[vwork->size[0] - 1] = vwork->size[0];
      }

      dim = 2;
      while (dim < i) {
        i2 = dim << 1;
        b_j = 1;
        for (pEnd = dim + 1; pEnd < i + 1; pEnd = qEnd + dim) {
          p = b_j;
          q = pEnd;
          qEnd = b_j + i2;
          if (qEnd > i + 1) {
            qEnd = i + 1;
          }

          k = 0;
          kEnd = qEnd - b_j;
          while (k + 1 <= kEnd) {
            if (sortLE(vwork, iidx->data[p - 1], iidx->data[q - 1])) {
              iwork->data[k] = iidx->data[p - 1];
              p++;
              if (p == pEnd) {
                while (q < qEnd) {
                  k++;
                  iwork->data[k] = iidx->data[q - 1];
                  q++;
                }
              }
            } else {
              iwork->data[k] = iidx->data[q - 1];
              q++;
              if (q == qEnd) {
                while (p < pEnd) {
                  k++;
                  iwork->data[k] = iidx->data[p - 1];
                  p++;
                }
              }
            }

            k++;
          }

          for (k = 0; k < kEnd; k++) {
            iidx->data[(b_j + k) - 1] = iwork->data[k];
          }

          b_j = qEnd;
        }

        dim = i2;
      }

      i = xwork->size[0];
      xwork->size[0] = vwork->size[0];
      emxEnsureCapacity_creal_T(xwork, i);
      for (k = 0; k <= n - 2; k++) {
        xwork->data[k] = vwork->data[k];
      }

      for (k = 0; k <= n - 2; k++) {
        vwork->data[k] = xwork->data[iidx->data[k] - 1];
      }
    }

    for (k = 0; k <= vlen; k++) {
      i = j + k * vstride;
      x->data[i] = vwork->data[k];
      idx->data[i] = iidx->data[k];
    }
  }

  emxFree_creal_T(&xwork);
  emxFree_int32_T(&iwork);
  emxFree_int32_T(&iidx);
  emxFree_creal_T(&vwork);
}

//
// Arguments    : emxArray_real_T *x
// Return Type  : void
//
void sort(emxArray_real_T *x)
{
  int dim;
  emxArray_real_T *vwork;
  int j;
  int vlen;
  int vstride;
  int k;
  emxArray_int32_T *b_vwork;
  dim = 0;
  if (x->size[0] != 1) {
    dim = -1;
  }

  emxInit_real_T(&vwork, 1);
  if (dim + 2 <= 1) {
    j = x->size[0];
  } else {
    j = 1;
  }

  vlen = j - 1;
  vstride = vwork->size[0];
  vwork->size[0] = j;
  emxEnsureCapacity_real_T(vwork, vstride);
  vstride = 1;
  for (k = 0; k <= dim; k++) {
    vstride *= x->size[0];
  }

  emxInit_int32_T(&b_vwork, 1);
  for (j = 0; j < vstride; j++) {
    for (k = 0; k <= vlen; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    sortIdx(vwork, b_vwork);
    for (k = 0; k <= vlen; k++) {
      x->data[j + k * vstride] = vwork->data[k];
    }
  }

  emxFree_int32_T(&b_vwork);
  emxFree_real_T(&vwork);
}

//
// File trailer for sort.cpp
//
// [EOF]
//
