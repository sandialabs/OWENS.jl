//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: insertionsort.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 15:21:39
//

// Include Files
#include "insertionsort.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// Arguments    : emxArray_int32_T *x
//                int xstart
//                int xend
// Return Type  : void
//
void insertionsort(emxArray_int32_T *x, int xstart, int xend)
{
  int i;
  int k;
  int xc;
  int idx;
  i = xstart + 1;
  for (k = i; k <= xend; k++) {
    xc = x->data[k - 1];
    idx = k - 1;
    while ((idx >= xstart) && (xc < x->data[idx - 1])) {
      x->data[idx] = x->data[idx - 1];
      idx--;
    }

    x->data[idx] = xc;
  }
}

//
// File trailer for insertionsort.cpp
//
// [EOF]
//
