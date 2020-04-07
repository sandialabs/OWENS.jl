//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: mldivide.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "mldivide.h"
#include "qrsolve.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include "xgeqp3.h"
#include "xgetrf.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_real_T *A
//                emxArray_real_T *B
// Return Type  : void
//
void b_mldivide(const emxArray_real_T *A, emxArray_real_T *B)
{
  emxArray_real_T *b_A;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  emxArray_real_T *b_B;
  int i;
  int mn;
  int loop_ub;
  int rankA;
  int m;
  double wj;
  int b_i;
  emxInit_real_T(&b_A, 2);
  emxInit_real_T(&tau, 1);
  emxInit_int32_T(&jpvt, 2);
  emxInit_real_T(&b_B, 1);
  if ((A->size[0] == 0) || (A->size[1] == 0) || (B->size[0] == 0)) {
    i = B->size[0];
    B->size[0] = A->size[1];
    emxEnsureCapacity_real_T(B, i);
    loop_ub = A->size[1];
    for (i = 0; i < loop_ub; i++) {
      B->data[i] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    mn = A->size[1];
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, i);
    loop_ub = A->size[0] * A->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_A->data[i] = A->data[i];
    }

    xgetrf(A->size[1], A->size[1], b_A, A->size[1], jpvt, &loop_ub);
    i = A->size[1];
    for (loop_ub = 0; loop_ub <= i - 2; loop_ub++) {
      m = jpvt->data[loop_ub];
      if (m != loop_ub + 1) {
        wj = B->data[loop_ub];
        B->data[loop_ub] = B->data[m - 1];
        B->data[m - 1] = wj;
      }
    }

    for (loop_ub = 0; loop_ub < mn; loop_ub++) {
      m = mn * loop_ub;
      if (B->data[loop_ub] != 0.0) {
        i = loop_ub + 2;
        for (b_i = i; b_i <= mn; b_i++) {
          B->data[b_i - 1] -= B->data[loop_ub] * b_A->data[(b_i + m) - 1];
        }
      }
    }

    for (loop_ub = mn; loop_ub >= 1; loop_ub--) {
      m = mn * (loop_ub - 1);
      wj = B->data[loop_ub - 1];
      if (wj != 0.0) {
        B->data[loop_ub - 1] = wj / b_A->data[(loop_ub + m) - 1];
        for (b_i = 0; b_i <= loop_ub - 2; b_i++) {
          B->data[b_i] -= B->data[loop_ub - 1] * b_A->data[b_i + m];
        }
      }
    }
  } else {
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, i);
    loop_ub = A->size[0] * A->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_A->data[i] = A->data[i];
    }

    xgeqp3(b_A, tau, jpvt);
    rankA = rankFromQR(b_A);
    i = b_B->size[0];
    b_B->size[0] = B->size[0];
    emxEnsureCapacity_real_T(b_B, i);
    loop_ub = B->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_B->data[i] = B->data[i];
    }

    i = B->size[0];
    B->size[0] = b_A->size[1];
    emxEnsureCapacity_real_T(B, i);
    loop_ub = b_A->size[1];
    for (i = 0; i < loop_ub; i++) {
      B->data[i] = 0.0;
    }

    m = b_A->size[0];
    loop_ub = b_A->size[0];
    mn = b_A->size[1];
    if (loop_ub < mn) {
      mn = loop_ub;
    }

    for (loop_ub = 0; loop_ub < mn; loop_ub++) {
      if (tau->data[loop_ub] != 0.0) {
        wj = b_B->data[loop_ub];
        i = loop_ub + 2;
        for (b_i = i; b_i <= m; b_i++) {
          wj += b_A->data[(b_i + b_A->size[0] * loop_ub) - 1] * b_B->data[b_i -
            1];
        }

        wj *= tau->data[loop_ub];
        if (wj != 0.0) {
          b_B->data[loop_ub] -= wj;
          i = loop_ub + 2;
          for (b_i = i; b_i <= m; b_i++) {
            b_B->data[b_i - 1] -= b_A->data[(b_i + b_A->size[0] * loop_ub) - 1] *
              wj;
          }
        }
      }
    }

    for (b_i = 0; b_i < rankA; b_i++) {
      B->data[jpvt->data[b_i] - 1] = b_B->data[b_i];
    }

    for (loop_ub = rankA; loop_ub >= 1; loop_ub--) {
      i = jpvt->data[loop_ub - 1];
      B->data[i - 1] /= b_A->data[(loop_ub + b_A->size[0] * (loop_ub - 1)) - 1];
      for (b_i = 0; b_i <= loop_ub - 2; b_i++) {
        B->data[jpvt->data[b_i] - 1] -= B->data[jpvt->data[loop_ub - 1] - 1] *
          b_A->data[b_i + b_A->size[0] * (loop_ub - 1)];
      }
    }
  }

  emxFree_real_T(&b_B);
  emxFree_int32_T(&jpvt);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
}

//
// Arguments    : const emxArray_real_T *A
//                const emxArray_real_T *B
//                emxArray_real_T *Y
// Return Type  : void
//
void mldivide(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *Y)
{
  emxArray_real_T *b_A;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  emxArray_real_T *b_B;
  int i;
  int n;
  int jBcol;
  int nb;
  int j;
  int mn;
  int m;
  int k;
  double wj;
  int b_nb;
  int b_i;
  emxInit_real_T(&b_A, 2);
  emxInit_real_T(&tau, 1);
  emxInit_int32_T(&jpvt, 2);
  emxInit_real_T(&b_B, 2);
  if ((A->size[0] == 0) || (A->size[1] == 0) || ((B->size[0] == 0) || (B->size[1]
        == 0))) {
    i = Y->size[0] * Y->size[1];
    Y->size[0] = A->size[1];
    Y->size[1] = B->size[1];
    emxEnsureCapacity_real_T(Y, i);
    jBcol = A->size[1] * B->size[1];
    for (i = 0; i < jBcol; i++) {
      Y->data[i] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    n = A->size[1];
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, i);
    jBcol = A->size[0] * A->size[1];
    for (i = 0; i < jBcol; i++) {
      b_A->data[i] = A->data[i];
    }

    xgetrf(A->size[1], A->size[1], b_A, A->size[1], jpvt, &jBcol);
    nb = B->size[1] - 1;
    i = Y->size[0] * Y->size[1];
    Y->size[0] = B->size[0];
    Y->size[1] = B->size[1];
    emxEnsureCapacity_real_T(Y, i);
    jBcol = B->size[0] * B->size[1];
    for (i = 0; i < jBcol; i++) {
      Y->data[i] = B->data[i];
    }

    i = A->size[1];
    for (jBcol = 0; jBcol <= i - 2; jBcol++) {
      mn = jpvt->data[jBcol];
      if (mn != jBcol + 1) {
        for (m = 0; m <= nb; m++) {
          wj = Y->data[jBcol + Y->size[0] * m];
          Y->data[jBcol + Y->size[0] * m] = Y->data[(mn + Y->size[0] * m) - 1];
          Y->data[(mn + Y->size[0] * m) - 1] = wj;
        }
      }
    }

    for (j = 0; j <= nb; j++) {
      jBcol = n * j;
      for (k = 0; k < n; k++) {
        m = n * k;
        i = k + jBcol;
        if (Y->data[i] != 0.0) {
          mn = k + 2;
          for (b_i = mn; b_i <= n; b_i++) {
            b_nb = (b_i + jBcol) - 1;
            Y->data[b_nb] -= Y->data[i] * b_A->data[(b_i + m) - 1];
          }
        }
      }
    }

    for (j = 0; j <= nb; j++) {
      jBcol = n * j - 1;
      for (k = n; k >= 1; k--) {
        m = n * (k - 1) - 1;
        i = k + jBcol;
        if (Y->data[i] != 0.0) {
          Y->data[i] /= b_A->data[k + m];
          for (b_i = 0; b_i <= k - 2; b_i++) {
            mn = (b_i + jBcol) + 1;
            Y->data[mn] -= Y->data[i] * b_A->data[(b_i + m) + 1];
          }
        }
      }
    }
  } else {
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, i);
    jBcol = A->size[0] * A->size[1];
    for (i = 0; i < jBcol; i++) {
      b_A->data[i] = A->data[i];
    }

    xgeqp3(b_A, tau, jpvt);
    n = rankFromQR(b_A);
    i = b_B->size[0] * b_B->size[1];
    b_B->size[0] = B->size[0];
    b_B->size[1] = B->size[1];
    emxEnsureCapacity_real_T(b_B, i);
    jBcol = B->size[0] * B->size[1];
    for (i = 0; i < jBcol; i++) {
      b_B->data[i] = B->data[i];
    }

    nb = B->size[1];
    i = Y->size[0] * Y->size[1];
    Y->size[0] = b_A->size[1];
    Y->size[1] = B->size[1];
    emxEnsureCapacity_real_T(Y, i);
    jBcol = b_A->size[1] * B->size[1];
    for (i = 0; i < jBcol; i++) {
      Y->data[i] = 0.0;
    }

    m = b_A->size[0];
    b_nb = B->size[1];
    jBcol = b_A->size[0];
    mn = b_A->size[1];
    if (jBcol < mn) {
      mn = jBcol;
    }

    for (j = 0; j < mn; j++) {
      if (tau->data[j] != 0.0) {
        for (k = 0; k < b_nb; k++) {
          wj = b_B->data[j + b_B->size[0] * k];
          i = j + 2;
          for (b_i = i; b_i <= m; b_i++) {
            wj += b_A->data[(b_i + b_A->size[0] * j) - 1] * b_B->data[(b_i +
              b_B->size[0] * k) - 1];
          }

          wj *= tau->data[j];
          if (wj != 0.0) {
            b_B->data[j + b_B->size[0] * k] -= wj;
            i = j + 2;
            for (b_i = i; b_i <= m; b_i++) {
              b_B->data[(b_i + b_B->size[0] * k) - 1] -= b_A->data[(b_i +
                b_A->size[0] * j) - 1] * wj;
            }
          }
        }
      }
    }

    for (k = 0; k < nb; k++) {
      for (b_i = 0; b_i < n; b_i++) {
        Y->data[(jpvt->data[b_i] + Y->size[0] * k) - 1] = b_B->data[b_i +
          b_B->size[0] * k];
      }

      for (j = n; j >= 1; j--) {
        i = jpvt->data[j - 1];
        Y->data[(i + Y->size[0] * k) - 1] /= b_A->data[(j + b_A->size[0] * (j -
          1)) - 1];
        for (b_i = 0; b_i <= j - 2; b_i++) {
          Y->data[(jpvt->data[b_i] + Y->size[0] * k) - 1] -= Y->data[(jpvt->
            data[j - 1] + Y->size[0] * k) - 1] * b_A->data[b_i + b_A->size[0] *
            (j - 1)];
        }
      }
    }
  }

  emxFree_real_T(&b_B);
  emxFree_int32_T(&jpvt);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
}

//
// File trailer for mldivide.cpp
//
// [EOF]
//
