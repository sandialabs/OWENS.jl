//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: xgetrf.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "xgetrf.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cmath>
#include <string.h>

// Function Definitions

//
// Arguments    : int m
//                int n
//                emxArray_real_T *A
//                int lda
//                emxArray_int32_T *ipiv
//                int *info
// Return Type  : void
//
void xgetrf(int m, int n, emxArray_real_T *A, int lda, emxArray_int32_T *ipiv,
            int *info)
{
  int yk;
  int b_n;
  int i;
  int jA;
  int u0;
  int j;
  int mmj;
  int b_tmp;
  int jp1j;
  int ix;
  double smax;
  double s;
  int i1;
  int ijA;
  if (m < n) {
    yk = m;
  } else {
    yk = n;
  }

  if (yk < 1) {
    b_n = 0;
  } else {
    b_n = yk;
  }

  i = ipiv->size[0] * ipiv->size[1];
  ipiv->size[0] = 1;
  ipiv->size[1] = b_n;
  emxEnsureCapacity_int32_T(ipiv, i);
  if (b_n > 0) {
    ipiv->data[0] = 1;
    yk = 1;
    for (jA = 2; jA <= b_n; jA++) {
      yk++;
      ipiv->data[jA - 1] = yk;
    }
  }

  *info = 0;
  if ((m >= 1) && (n >= 1)) {
    u0 = m - 1;
    if (u0 >= n) {
      u0 = n;
    }

    for (j = 0; j < u0; j++) {
      mmj = m - j;
      b_tmp = j * (lda + 1);
      jp1j = b_tmp + 2;
      if (mmj < 1) {
        yk = -1;
      } else {
        yk = 0;
        if (mmj > 1) {
          ix = b_tmp;
          smax = std::abs(A->data[b_tmp]);
          for (jA = 2; jA <= mmj; jA++) {
            ix++;
            s = std::abs(A->data[ix]);
            if (s > smax) {
              yk = jA - 1;
              smax = s;
            }
          }
        }
      }

      if (A->data[b_tmp + yk] != 0.0) {
        if (yk != 0) {
          yk += j;
          ipiv->data[j] = yk + 1;
          ix = j;
          for (jA = 0; jA < n; jA++) {
            smax = A->data[ix];
            A->data[ix] = A->data[yk];
            A->data[yk] = smax;
            ix += lda;
            yk += lda;
          }
        }

        i = b_tmp + mmj;
        for (yk = jp1j; yk <= i; yk++) {
          A->data[yk - 1] /= A->data[b_tmp];
        }
      } else {
        *info = j + 1;
      }

      b_n = n - j;
      yk = b_tmp + lda;
      jA = yk;
      for (jp1j = 0; jp1j <= b_n - 2; jp1j++) {
        smax = A->data[yk];
        if (A->data[yk] != 0.0) {
          ix = b_tmp + 1;
          i = jA + 2;
          i1 = mmj + jA;
          for (ijA = i; ijA <= i1; ijA++) {
            A->data[ijA - 1] += A->data[ix] * -smax;
            ix++;
          }
        }

        yk += lda;
        jA += lda;
      }
    }

    if ((*info == 0) && (m <= n) && (!(A->data[(m + A->size[0] * (m - 1)) - 1]
          != 0.0))) {
      *info = m;
    }
  }
}

//
// File trailer for xgetrf.cpp
//
// [EOF]
//
