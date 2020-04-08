//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: schur.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:47:29
//

// Include Files
#include "schur.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include "xdhseqr.h"
#include "xzlarf.h"
#include "xzlarfg.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_real_T *A
//                emxArray_real_T *V
//                emxArray_real_T *T
// Return Type  : void
//
void schur(const emxArray_real_T *A, emxArray_real_T *V, emxArray_real_T *T)
{
  int nx;
  boolean_T p;
  int istart;
  int n;
  int i;
  emxArray_real_T *tau;
  emxArray_real_T *work;
  int b_n;
  int iaii;
  int j;
  int b_i;
  int im1n;
  int in;
  double alpha1;
  double temp;
  int nh;
  int jy;
  int ia;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int ix;
  int i1;
  int exitg1;
  nx = A->size[0] * A->size[1];
  p = true;
  for (istart = 0; istart < nx; istart++) {
    if ((!p) || (rtIsInf(A->data[istart]) || rtIsNaN(A->data[istart]))) {
      p = false;
    }
  }

  if (!p) {
    i = V->size[0] * V->size[1];
    V->size[0] = A->size[0];
    V->size[1] = A->size[1];
    emxEnsureCapacity_real_T(V, i);
    nx = A->size[0] * A->size[1];
    for (i = 0; i < nx; i++) {
      V->data[i] = rtNaN;
    }

    nx = V->size[0];
    if ((V->size[0] != 0) && (V->size[1] != 0) && (1 < V->size[0])) {
      istart = 2;
      if (V->size[0] - 2 < V->size[1] - 1) {
        iaii = V->size[0] - 1;
      } else {
        iaii = V->size[1];
      }

      for (j = 0; j < iaii; j++) {
        for (b_i = istart; b_i <= nx; b_i++) {
          V->data[(b_i + V->size[0] * j) - 1] = 0.0;
        }

        istart++;
      }
    }

    i = T->size[0] * T->size[1];
    T->size[0] = A->size[0];
    T->size[1] = A->size[1];
    emxEnsureCapacity_real_T(T, i);
    nx = A->size[0] * A->size[1];
    for (i = 0; i < nx; i++) {
      T->data[i] = rtNaN;
    }
  } else {
    n = A->size[0];
    i = T->size[0] * T->size[1];
    T->size[0] = A->size[0];
    T->size[1] = A->size[1];
    emxEnsureCapacity_real_T(T, i);
    nx = A->size[0] * A->size[1];
    for (i = 0; i < nx; i++) {
      T->data[i] = A->data[i];
    }

    emxInit_real_T(&tau, 1);
    emxInit_real_T(&work, 1);
    b_n = A->size[0];
    i = tau->size[0];
    if (A->size[0] < 1) {
      tau->size[0] = 0;
    } else {
      tau->size[0] = A->size[0] - 1;
    }

    emxEnsureCapacity_real_T(tau, i);
    i = work->size[0];
    work->size[0] = A->size[0];
    emxEnsureCapacity_real_T(work, i);
    nx = A->size[0];
    for (i = 0; i < nx; i++) {
      work->data[i] = 0.0;
    }

    i = A->size[0];
    for (b_i = 0; b_i <= i - 2; b_i++) {
      im1n = b_i * b_n + 2;
      in = (b_i + 1) * b_n;
      alpha1 = T->data[(b_i + T->size[0] * b_i) + 1];
      nx = b_i + 3;
      if (nx >= b_n) {
        nx = b_n;
      }

      temp = xzlarfg((b_n - b_i) - 1, &alpha1, T, nx + b_i * b_n);
      tau->data[b_i] = temp;
      T->data[(b_i + T->size[0] * b_i) + 1] = 1.0;
      nx = (b_n - b_i) - 3;
      jy = (b_i + im1n) - 1;
      iaii = in + 1;
      if (temp != 0.0) {
        lastv = nx + 1;
        nx += jy;
        while ((lastv + 1 > 0) && (T->data[nx + 1] == 0.0)) {
          lastv--;
          nx--;
        }

        lastc = b_n;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          nx = in + lastc;
          ia = nx;
          do {
            exitg1 = 0;
            if ((b_n > 0) && (ia <= nx + lastv * b_n)) {
              if (T->data[ia - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ia += b_n;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = -1;
        lastc = 0;
      }

      if (lastv + 1 > 0) {
        if (lastc != 0) {
          for (nx = 0; nx < lastc; nx++) {
            work->data[nx] = 0.0;
          }

          ix = jy;
          i1 = (in + b_n * lastv) + 1;
          for (istart = iaii; b_n < 0 ? istart >= i1 : istart <= i1; istart +=
               b_n) {
            nx = 0;
            nh = (istart + lastc) - 1;
            for (ia = istart; ia <= nh; ia++) {
              work->data[nx] += T->data[ia - 1] * T->data[ix];
              nx++;
            }

            ix++;
          }
        }

        if (!(-tau->data[b_i] == 0.0)) {
          nx = in;
          for (j = 0; j <= lastv; j++) {
            if (T->data[jy] != 0.0) {
              temp = T->data[jy] * -tau->data[b_i];
              ix = 0;
              i1 = nx + 1;
              nh = lastc + nx;
              for (istart = i1; istart <= nh; istart++) {
                T->data[istart - 1] += work->data[ix] * temp;
                ix++;
              }
            }

            jy++;
            nx += b_n;
          }
        }
      }

      xzlarf((b_n - b_i) - 1, (b_n - b_i) - 1, b_i + im1n, tau->data[b_i], T,
             (b_i + in) + 2, b_n, work);
      T->data[(b_i + T->size[0] * b_i) + 1] = alpha1;
    }

    i = V->size[0] * V->size[1];
    V->size[0] = T->size[0];
    V->size[1] = T->size[1];
    emxEnsureCapacity_real_T(V, i);
    nx = T->size[0] * T->size[1];
    for (i = 0; i < nx; i++) {
      V->data[i] = T->data[i];
    }

    if (A->size[0] != 0) {
      nh = A->size[0] - 1;
      for (j = n; j >= 2; j--) {
        ia = (j - 1) * n - 1;
        for (b_i = 0; b_i <= j - 2; b_i++) {
          V->data[(ia + b_i) + 1] = 0.0;
        }

        nx = ia - n;
        i = j + 1;
        for (b_i = i; b_i <= n; b_i++) {
          V->data[ia + b_i] = V->data[nx + b_i];
        }

        i = n + 1;
        for (b_i = i; b_i <= n; b_i++) {
          V->data[ia + b_i] = 0.0;
        }
      }

      for (b_i = 0; b_i < n; b_i++) {
        V->data[b_i] = 0.0;
      }

      V->data[0] = 1.0;
      i = A->size[0] + 1;
      for (j = i; j <= n; j++) {
        ia = (j - 1) * n;
        for (b_i = 0; b_i < n; b_i++) {
          V->data[ia + b_i] = 0.0;
        }

        V->data[(ia + j) - 1] = 1.0;
      }

      if (A->size[0] - 1 >= 1) {
        i = A->size[0] - 2;
        for (j = nh; j <= i; j++) {
          ia = (n + j * n) + 1;
          i1 = n - 2;
          for (b_i = 0; b_i <= i1; b_i++) {
            V->data[ia + b_i] = 0.0;
          }

          V->data[ia + j] = 1.0;
        }

        ix = A->size[0] - 2;
        i = work->size[0];
        work->size[0] = V->size[1];
        emxEnsureCapacity_real_T(work, i);
        nx = V->size[1];
        for (i = 0; i < nx; i++) {
          work->data[i] = 0.0;
        }

        for (b_i = A->size[0] - 1; b_i >= 1; b_i--) {
          iaii = (n + b_i) + (b_i - 1) * n;
          if (b_i < n - 1) {
            V->data[iaii] = 1.0;
            xzlarf(n - b_i, nh - b_i, iaii + 1, tau->data[ix], V, (iaii + n) + 1,
                   n, work);
            nx = iaii + 2;
            i = (iaii + n) - b_i;
            for (istart = nx; istart <= i; istart++) {
              V->data[istart - 1] *= -tau->data[ix];
            }
          }

          V->data[iaii] = 1.0 - tau->data[ix];
          for (j = 0; j <= b_i - 2; j++) {
            V->data[(iaii - j) - 1] = 0.0;
          }

          ix--;
        }
      }
    }

    emxFree_real_T(&work);
    emxFree_real_T(&tau);
    eml_dlahqr(T, V);
    nx = T->size[0];
    if ((T->size[0] != 0) && (T->size[1] != 0) && (3 < T->size[0])) {
      istart = 4;
      if (T->size[0] - 4 < T->size[1] - 1) {
        iaii = T->size[0] - 3;
      } else {
        iaii = T->size[1];
      }

      for (j = 0; j < iaii; j++) {
        for (b_i = istart; b_i <= nx; b_i++) {
          T->data[(b_i + T->size[0] * j) - 1] = 0.0;
        }

        istart++;
      }
    }
  }
}

//
// File trailer for schur.cpp
//
// [EOF]
//
