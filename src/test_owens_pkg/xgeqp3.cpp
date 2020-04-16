//
// File: xgeqp3.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//

// Include Files
#include "xgeqp3.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include "xgerc.h"
#include "xnrm2.h"
#include "xzlarfg.h"
#include <cmath>
#include <string.h>

// Function Definitions

//
// Arguments    : emxArray_real_T *A
//                emxArray_real_T *tau
//                emxArray_int32_T *jpvt
// Return Type  : void
//
void xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, emxArray_int32_T *jpvt)
{
  int n;
  int itemp;
  int minmana;
  int i;
  emxArray_real_T *work;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  boolean_T guard1 = false;
  int m;
  int ma;
  int minmn;
  int b_i;
  double d;
  int ip1;
  int iy;
  int ii;
  int nmi;
  int mmi;
  int ix;
  double smax;
  double s;
  int ic0;
  int lastv;
  boolean_T exitg2;
  int ia;
  int exitg1;
  n = A->size[1] - 1;
  itemp = A->size[0];
  minmana = A->size[1];
  if (itemp < minmana) {
    minmana = itemp;
  }

  i = tau->size[0];
  tau->size[0] = minmana;
  emxEnsureCapacity_real_T(tau, i);
  for (i = 0; i < minmana; i++) {
    tau->data[i] = 0.0;
  }

  emxInit_real_T(&work, 1);
  emxInit_real_T(&vn1, 1);
  emxInit_real_T(&vn2, 1);
  guard1 = false;
  if ((A->size[0] == 0) || (A->size[1] == 0)) {
    guard1 = true;
  } else {
    itemp = A->size[0];
    minmana = A->size[1];
    if (itemp < minmana) {
      minmana = itemp;
    }

    if (minmana < 1) {
      guard1 = true;
    } else {
      i = jpvt->size[0] * jpvt->size[1];
      jpvt->size[0] = 1;
      jpvt->size[1] = A->size[1];
      emxEnsureCapacity_int32_T(jpvt, i);
      itemp = A->size[1];
      for (i = 0; i < itemp; i++) {
        jpvt->data[i] = 0;
      }

      for (itemp = 0; itemp <= n; itemp++) {
        jpvt->data[itemp] = itemp + 1;
      }

      m = A->size[0];
      n = A->size[1];
      ma = A->size[0];
      if (m < n) {
        minmn = m;
      } else {
        minmn = n;
      }

      i = work->size[0];
      work->size[0] = A->size[1];
      emxEnsureCapacity_real_T(work, i);
      itemp = A->size[1];
      for (i = 0; i < itemp; i++) {
        work->data[i] = 0.0;
      }

      i = vn1->size[0];
      vn1->size[0] = A->size[1];
      emxEnsureCapacity_real_T(vn1, i);
      itemp = A->size[1];
      for (i = 0; i < itemp; i++) {
        vn1->data[i] = 0.0;
      }

      i = vn2->size[0];
      vn2->size[0] = A->size[1];
      emxEnsureCapacity_real_T(vn2, i);
      itemp = A->size[1];
      for (i = 0; i < itemp; i++) {
        vn2->data[i] = 0.0;
      }

      for (itemp = 0; itemp < n; itemp++) {
        d = xnrm2(m, A, itemp * ma + 1);
        vn1->data[itemp] = d;
        vn2->data[itemp] = d;
      }

      for (b_i = 0; b_i < minmn; b_i++) {
        ip1 = b_i + 2;
        iy = b_i * ma;
        ii = iy + b_i;
        nmi = n - b_i;
        mmi = (m - b_i) - 1;
        if (nmi < 1) {
          minmana = -1;
        } else {
          minmana = 0;
          if (nmi > 1) {
            ix = b_i;
            smax = std::abs(vn1->data[b_i]);
            for (itemp = 2; itemp <= nmi; itemp++) {
              ix++;
              s = std::abs(vn1->data[ix]);
              if (s > smax) {
                minmana = itemp - 1;
                smax = s;
              }
            }
          }
        }

        minmana += b_i;
        if (minmana + 1 != b_i + 1) {
          ix = minmana * ma;
          for (itemp = 0; itemp < m; itemp++) {
            smax = A->data[ix];
            A->data[ix] = A->data[iy];
            A->data[iy] = smax;
            ix++;
            iy++;
          }

          itemp = jpvt->data[minmana];
          jpvt->data[minmana] = jpvt->data[b_i];
          jpvt->data[b_i] = itemp;
          vn1->data[minmana] = vn1->data[b_i];
          vn2->data[minmana] = vn2->data[b_i];
        }

        if (b_i + 1 < m) {
          s = A->data[ii];
          d = xzlarfg(mmi + 1, &s, A, ii + 2);
          tau->data[b_i] = d;
          A->data[ii] = s;
        } else {
          d = 0.0;
          tau->data[b_i] = 0.0;
        }

        if (b_i + 1 < n) {
          s = A->data[ii];
          A->data[ii] = 1.0;
          ic0 = (ii + ma) + 1;
          if (d != 0.0) {
            lastv = mmi + 1;
            minmana = ii + mmi;
            while ((lastv > 0) && (A->data[minmana] == 0.0)) {
              lastv--;
              minmana--;
            }

            nmi--;
            exitg2 = false;
            while ((!exitg2) && (nmi > 0)) {
              minmana = ic0 + (nmi - 1) * ma;
              ia = minmana;
              do {
                exitg1 = 0;
                if (ia <= (minmana + lastv) - 1) {
                  if (A->data[ia - 1] != 0.0) {
                    exitg1 = 1;
                  } else {
                    ia++;
                  }
                } else {
                  nmi--;
                  exitg1 = 2;
                }
              } while (exitg1 == 0);

              if (exitg1 == 1) {
                exitg2 = true;
              }
            }
          } else {
            lastv = 0;
            nmi = 0;
          }

          if (lastv > 0) {
            if (nmi != 0) {
              for (iy = 0; iy < nmi; iy++) {
                work->data[iy] = 0.0;
              }

              iy = 0;
              i = ic0 + ma * (nmi - 1);
              for (minmana = ic0; ma < 0 ? minmana >= i : minmana <= i; minmana +=
                   ma) {
                ix = ii;
                smax = 0.0;
                itemp = (minmana + lastv) - 1;
                for (ia = minmana; ia <= itemp; ia++) {
                  smax += A->data[ia - 1] * A->data[ix];
                  ix++;
                }

                work->data[iy] += smax;
                iy++;
              }
            }

            xgerc(lastv, nmi, -tau->data[b_i], ii + 1, work, A, ic0, ma);
          }

          A->data[ii] = s;
        }

        for (itemp = ip1; itemp <= n; itemp++) {
          minmana = b_i + (itemp - 1) * ma;
          d = vn1->data[itemp - 1];
          if (d != 0.0) {
            smax = std::abs(A->data[minmana]) / d;
            smax = 1.0 - smax * smax;
            if (smax < 0.0) {
              smax = 0.0;
            }

            s = d / vn2->data[itemp - 1];
            s = smax * (s * s);
            if (s <= 1.4901161193847656E-8) {
              if (b_i + 1 < m) {
                d = xnrm2(mmi, A, minmana + 2);
                vn1->data[itemp - 1] = d;
                vn2->data[itemp - 1] = d;
              } else {
                vn1->data[itemp - 1] = 0.0;
                vn2->data[itemp - 1] = 0.0;
              }
            } else {
              vn1->data[itemp - 1] = d * std::sqrt(smax);
            }
          }
        }
      }
    }
  }

  if (guard1) {
    i = jpvt->size[0] * jpvt->size[1];
    jpvt->size[0] = 1;
    jpvt->size[1] = A->size[1];
    emxEnsureCapacity_int32_T(jpvt, i);
    itemp = A->size[1];
    for (i = 0; i < itemp; i++) {
      jpvt->data[i] = 0;
    }

    for (itemp = 0; itemp <= n; itemp++) {
      jpvt->data[itemp] = itemp + 1;
    }
  }

  emxFree_real_T(&vn2);
  emxFree_real_T(&vn1);
  emxFree_real_T(&work);
}

//
// File trailer for xgeqp3.cpp
//
// [EOF]
//
