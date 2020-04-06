//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: xzlarf.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//

// Include Files
#include "xzlarf.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "xgerc.h"
#include <string.h>

// Function Definitions

//
// Arguments    : int m
//                int n
//                int iv0
//                double tau
//                emxArray_real_T *C
//                int ic0
//                int ldc
//                emxArray_real_T *work
// Return Type  : void
//
void xzlarf(int m, int n, int iv0, double tau, emxArray_real_T *C, int ic0, int
            ldc, emxArray_real_T *work)
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  int b_i;
  int iac;
  int ia;
  int ix;
  int exitg1;
  double c;
  int i1;
  if (tau != 0.0) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C->data[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      i = ic0 + (lastc - 1) * ldc;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C->data[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
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
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      for (i = 0; i < lastc; i++) {
        work->data[i] = 0.0;
      }

      i = 0;
      b_i = ic0 + ldc * (lastc - 1);
      for (iac = ic0; ldc < 0 ? iac >= b_i : iac <= b_i; iac += ldc) {
        ix = iv0;
        c = 0.0;
        i1 = (iac + lastv) - 1;
        for (ia = iac; ia <= i1; ia++) {
          c += C->data[ia - 1] * C->data[ix - 1];
          ix++;
        }

        work->data[i] += c;
        i++;
      }
    }

    xgerc(lastv, lastc, -tau, iv0, work, C, ic0, ldc);
  }
}

//
// File trailer for xzlarf.cpp
//
// [EOF]
//
