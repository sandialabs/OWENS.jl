//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: linear_interp.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "linear_interp.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_real_T *x_array
//                const emxArray_real_T *y_array
//                const emxArray_real_T *xnew
//                emxArray_real_T *ynew
// Return Type  : void
//
void linear_interp(const emxArray_real_T *x_array, const emxArray_real_T
                   *y_array, const emxArray_real_T *xnew, emxArray_real_T *ynew)
{
  unsigned int unnamed_idx_1;
  int i;
  int idx;
  int b_i;
  double b_xnew;
  int c_i;
  int n;
  double min_x;
  int k;
  boolean_T exitg1;
  double max_x;
  double d;
  unnamed_idx_1 = static_cast<unsigned int>(xnew->size[1]);
  i = ynew->size[0] * ynew->size[1];
  ynew->size[0] = 1;
  ynew->size[1] = static_cast<int>(unnamed_idx_1);
  emxEnsureCapacity_real_T(ynew, i);
  idx = static_cast<int>(unnamed_idx_1);
  for (i = 0; i < idx; i++) {
    ynew->data[i] = 0.0;
  }

  i = xnew->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    b_xnew = xnew->data[b_i];
    c_i = 0;
    n = x_array->size[1];
    if (x_array->size[1] <= 2) {
      if (x_array->size[1] == 1) {
        min_x = x_array->data[0];
      } else if ((x_array->data[0] > x_array->data[1]) || (rtIsNaN(x_array->
                   data[0]) && (!rtIsNaN(x_array->data[1])))) {
        min_x = x_array->data[1];
      } else {
        min_x = x_array->data[0];
      }
    } else {
      if (!rtIsNaN(x_array->data[0])) {
        idx = 1;
      } else {
        idx = 0;
        k = 2;
        exitg1 = false;
        while ((!exitg1) && (k <= x_array->size[1])) {
          if (!rtIsNaN(x_array->data[k - 1])) {
            idx = k;
            exitg1 = true;
          } else {
            k++;
          }
        }
      }

      if (idx == 0) {
        min_x = x_array->data[0];
      } else {
        min_x = x_array->data[idx - 1];
        idx++;
        for (k = idx; k <= n; k++) {
          d = x_array->data[k - 1];
          if (min_x > d) {
            min_x = d;
          }
        }
      }
    }

    n = x_array->size[1];
    if (x_array->size[1] <= 2) {
      if (x_array->size[1] == 1) {
        max_x = x_array->data[0];
      } else if ((x_array->data[0] < x_array->data[1]) || (rtIsNaN(x_array->
                   data[0]) && (!rtIsNaN(x_array->data[1])))) {
        max_x = x_array->data[1];
      } else {
        max_x = x_array->data[0];
      }
    } else {
      if (!rtIsNaN(x_array->data[0])) {
        idx = 1;
      } else {
        idx = 0;
        k = 2;
        exitg1 = false;
        while ((!exitg1) && (k <= x_array->size[1])) {
          if (!rtIsNaN(x_array->data[k - 1])) {
            idx = k;
            exitg1 = true;
          } else {
            k++;
          }
        }
      }

      if (idx == 0) {
        max_x = x_array->data[0];
      } else {
        max_x = x_array->data[idx - 1];
        idx++;
        for (k = idx; k <= n; k++) {
          d = x_array->data[k - 1];
          if (max_x < d) {
            max_x = d;
          }
        }
      }
    }

    d = xnew->data[b_i];
    if (!(d < min_x)) {
      if (d > max_x) {
        //      ynew = y_array(max_x_idx);
        c_i = x_array->size[1] - 2;
      } else {
        idx = 0;
        exitg1 = false;
        while ((!exitg1) && (idx <= x_array->size[1] - 2)) {
          if ((x_array->data[idx] <= b_xnew) && (x_array->data[idx + 1] >=
               b_xnew)) {
            c_i = idx;
            exitg1 = true;
          } else {
            idx++;
          }
        }
      }
    } else {
      //  extrapolate linearly at max and min
      //      ynew = y_array(min_x_idx);
    }

    ynew->data[b_i] = y_array->data[c_i] + (xnew->data[b_i] - x_array->data[c_i])
      / (x_array->data[c_i + 1] - x_array->data[c_i]) * (y_array->data[c_i + 1]
      - y_array->data[c_i]);
  }
}

//
// File trailer for linear_interp.cpp
//
// [EOF]
//
