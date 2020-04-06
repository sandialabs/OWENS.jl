//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: norm.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//

// Include Files
#include "norm.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cmath>
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_real_T *x
// Return Type  : double
//
double b_norm(const emxArray_real_T *x)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = 0.0;
    if (x->size[0] == 1) {
      y = std::abs(x->data[0]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = x->size[0];
      for (k = 0; k < kend; k++) {
        absxk = std::abs(x->data[k]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// File trailer for norm.cpp
//
// [EOF]
//
