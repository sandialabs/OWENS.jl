/*
 * File: mldivide.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "mldivide.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "test_transient1_emxutil.h"
#include "xzgeqp3.h"
#include <math.h>
#include <string.h>

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *A
 *                emxArray_real_T *B
 * Return Type  : void
 */
void mldiv(const emxArray_real_T *A, emxArray_real_T *B)
{
  emxArray_real_T *b_A;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  emxArray_real_T *b_B;
  int i;
  int n;
  int minmana;
  int jp1j;
  int maxmn;
  int rankR;
  int ldap1;
  int u1;
  int j;
  int mmj_tmp;
  int jj;
  int b_i;
  double tol;
  int ix;
  double s;
  emxInit_real_T(&b_A, 2);
  emxInit_real_T(&tau, 1);
  emxInit_int32_T(&jpvt, 2);
  emxInit_real_T(&b_B, 1);
  if ((A->size[0] == 0) || (A->size[1] == 0) || (B->size[0] == 0)) {
    i = B->size[0];
    B->size[0] = A->size[1];
    emxEnsureCapacity_real_T(B, i);
    minmana = A->size[1];
    for (i = 0; i < minmana; i++) {
      B->data[i] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    n = A->size[1];
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, i);
    minmana = A->size[0] * A->size[1];
    for (i = 0; i < minmana; i++) {
      b_A->data[i] = A->data[i];
    }

    minmana = A->size[1];
    i = jpvt->size[0] * jpvt->size[1];
    jpvt->size[0] = 1;
    jpvt->size[1] = A->size[1];
    emxEnsureCapacity_int32_T(jpvt, i);
    jpvt->data[0] = 1;
    maxmn = 1;
    for (rankR = 2; rankR <= minmana; rankR++) {
      maxmn++;
      jpvt->data[rankR - 1] = maxmn;
    }

    ldap1 = A->size[1];
    jp1j = A->size[1] - 1;
    u1 = A->size[1];
    if (jp1j < u1) {
      u1 = jp1j;
    }

    for (j = 0; j < u1; j++) {
      mmj_tmp = n - j;
      minmana = j * (n + 1);
      jj = j * (ldap1 + 1);
      jp1j = minmana + 2;
      if (mmj_tmp < 1) {
        maxmn = -1;
      } else {
        maxmn = 0;
        if (mmj_tmp > 1) {
          ix = minmana;
          tol = fabs(b_A->data[jj]);
          for (rankR = 2; rankR <= mmj_tmp; rankR++) {
            ix++;
            s = fabs(b_A->data[ix]);
            if (s > tol) {
              maxmn = rankR - 1;
              tol = s;
            }
          }
        }
      }

      if (b_A->data[jj + maxmn] != 0.0) {
        if (maxmn != 0) {
          maxmn += j;
          jpvt->data[j] = maxmn + 1;
          ix = j;
          for (rankR = 0; rankR < n; rankR++) {
            tol = b_A->data[ix];
            b_A->data[ix] = b_A->data[maxmn];
            b_A->data[maxmn] = tol;
            ix += n;
            maxmn += n;
          }
        }

        i = jj + mmj_tmp;
        for (b_i = jp1j; b_i <= i; b_i++) {
          b_A->data[b_i - 1] /= b_A->data[jj];
        }
      }

      maxmn = minmana + n;
      minmana = jj + ldap1;
      for (jp1j = 0; jp1j <= mmj_tmp - 2; jp1j++) {
        tol = b_A->data[maxmn];
        if (b_A->data[maxmn] != 0.0) {
          ix = jj + 1;
          i = minmana + 2;
          b_i = mmj_tmp + minmana;
          for (rankR = i; rankR <= b_i; rankR++) {
            b_A->data[rankR - 1] += b_A->data[ix] * -tol;
            ix++;
          }
        }

        maxmn += n;
        minmana += n;
      }
    }

    i = A->size[1];
    for (maxmn = 0; maxmn <= i - 2; maxmn++) {
      b_i = jpvt->data[maxmn];
      if (b_i != maxmn + 1) {
        tol = B->data[maxmn];
        B->data[maxmn] = B->data[b_i - 1];
        B->data[b_i - 1] = tol;
      }
    }

    for (rankR = 0; rankR < n; rankR++) {
      maxmn = n * rankR;
      if (B->data[rankR] != 0.0) {
        i = rankR + 2;
        for (b_i = i; b_i <= n; b_i++) {
          B->data[b_i - 1] -= B->data[rankR] * b_A->data[(b_i + maxmn) - 1];
        }
      }
    }

    for (rankR = n; rankR >= 1; rankR--) {
      maxmn = n * (rankR - 1);
      tol = B->data[rankR - 1];
      if (tol != 0.0) {
        B->data[rankR - 1] = tol / b_A->data[(rankR + maxmn) - 1];
        for (b_i = 0; b_i <= rankR - 2; b_i++) {
          B->data[b_i] -= B->data[rankR - 1] * b_A->data[b_i + maxmn];
        }
      }
    }
  } else {
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, i);
    minmana = A->size[0] * A->size[1];
    for (i = 0; i < minmana; i++) {
      b_A->data[i] = A->data[i];
    }

    n = A->size[1];
    jp1j = A->size[0];
    minmana = A->size[1];
    if (jp1j < minmana) {
      minmana = jp1j;
    }

    i = tau->size[0];
    tau->size[0] = minmana;
    emxEnsureCapacity_real_T(tau, i);
    for (i = 0; i < minmana; i++) {
      tau->data[i] = 0.0;
    }

    i = jpvt->size[0] * jpvt->size[1];
    jpvt->size[0] = 1;
    jpvt->size[1] = A->size[1];
    emxEnsureCapacity_int32_T(jpvt, i);
    minmana = A->size[1];
    for (i = 0; i < minmana; i++) {
      jpvt->data[i] = 0;
    }

    for (rankR = 0; rankR < n; rankR++) {
      jpvt->data[rankR] = rankR + 1;
    }

    qrpf(b_A, A->size[0], A->size[1], tau, jpvt);
    rankR = 0;
    if (b_A->size[0] < b_A->size[1]) {
      minmana = b_A->size[0];
      maxmn = b_A->size[1];
    } else {
      minmana = b_A->size[1];
      maxmn = b_A->size[0];
    }

    if (minmana > 0) {
      tol = fmin(1.4901161193847656E-8, 2.2204460492503131E-15 * (double)maxmn) *
        fabs(b_A->data[0]);
      while ((rankR < minmana) && (!(fabs(b_A->data[rankR + b_A->size[0] * rankR])
               <= tol))) {
        rankR++;
      }
    }

    i = b_B->size[0];
    b_B->size[0] = B->size[0];
    emxEnsureCapacity_real_T(b_B, i);
    minmana = B->size[0];
    for (i = 0; i < minmana; i++) {
      b_B->data[i] = B->data[i];
    }

    i = B->size[0];
    B->size[0] = b_A->size[1];
    emxEnsureCapacity_real_T(B, i);
    minmana = b_A->size[1];
    for (i = 0; i < minmana; i++) {
      B->data[i] = 0.0;
    }

    maxmn = b_A->size[0];
    jp1j = b_A->size[0];
    minmana = b_A->size[1];
    if (jp1j < minmana) {
      minmana = jp1j;
    }

    for (j = 0; j < minmana; j++) {
      if (tau->data[j] != 0.0) {
        tol = b_B->data[j];
        i = j + 2;
        for (b_i = i; b_i <= maxmn; b_i++) {
          tol += b_A->data[(b_i + b_A->size[0] * j) - 1] * b_B->data[b_i - 1];
        }

        tol *= tau->data[j];
        if (tol != 0.0) {
          b_B->data[j] -= tol;
          i = j + 2;
          for (b_i = i; b_i <= maxmn; b_i++) {
            b_B->data[b_i - 1] -= b_A->data[(b_i + b_A->size[0] * j) - 1] * tol;
          }
        }
      }
    }

    for (b_i = 0; b_i < rankR; b_i++) {
      B->data[jpvt->data[b_i] - 1] = b_B->data[b_i];
    }

    for (j = rankR; j >= 1; j--) {
      i = jpvt->data[j - 1];
      B->data[i - 1] /= b_A->data[(j + b_A->size[0] * (j - 1)) - 1];
      for (b_i = 0; b_i <= j - 2; b_i++) {
        B->data[jpvt->data[b_i] - 1] -= B->data[jpvt->data[j - 1] - 1] *
          b_A->data[b_i + b_A->size[0] * (j - 1)];
      }
    }
  }

  emxFree_real_T(&b_B);
  emxFree_int32_T(&jpvt);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
}

/*
 * File trailer for mldivide.c
 *
 * [EOF]
 */
