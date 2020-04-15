//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: qrsolve.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "qrsolve.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cmath>
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_real_T *A
// Return Type  : int
//
int rankFromQR(const emxArray_real_T *A)
{
  int b_r;
  int minmn;
  int maxmn;
  double tol;
  b_r = 0;
  if (A->size[0] < A->size[1]) {
    minmn = A->size[0];
    maxmn = A->size[1];
  } else {
    minmn = A->size[1];
    maxmn = A->size[0];
  }

  if (minmn > 0) {
    tol = 2.2204460492503131E-15 * static_cast<double>(maxmn);
    if (1.4901161193847656E-8 < tol) {
      tol = 1.4901161193847656E-8;
    }

    tol *= std::abs(A->data[0]);
    while ((b_r < minmn) && (!(std::abs(A->data[b_r + A->size[0] * b_r]) <= tol)))
    {
      b_r++;
    }
  }

  return b_r;
}

//
// File trailer for qrsolve.cpp
//
// [EOF]
//
