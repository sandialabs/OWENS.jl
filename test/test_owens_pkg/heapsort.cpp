//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: heapsort.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:47:29
//

// Include Files
#include "heapsort.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Declarations
static void heapify(emxArray_int32_T *x, int idx, int xstart, int xend);

// Function Definitions

//
// Arguments    : emxArray_int32_T *x
//                int idx
//                int xstart
//                int xend
// Return Type  : void
//
static void heapify(emxArray_int32_T *x, int idx, int xstart, int xend)
{
  boolean_T changed;
  int extremumIdx;
  int leftIdx;
  boolean_T exitg1;
  int extremum;
  int cmpIdx;
  int xcmp;
  changed = true;
  extremumIdx = (idx + xstart) - 2;
  leftIdx = ((idx << 1) + xstart) - 1;
  exitg1 = false;
  while ((!exitg1) && (leftIdx < xend)) {
    changed = false;
    extremum = x->data[extremumIdx];
    cmpIdx = leftIdx - 1;
    xcmp = x->data[leftIdx - 1];
    if (x->data[leftIdx - 1] < x->data[leftIdx]) {
      cmpIdx = leftIdx;
      xcmp = x->data[leftIdx];
    }

    if (x->data[extremumIdx] < xcmp) {
      x->data[extremumIdx] = xcmp;
      x->data[cmpIdx] = extremum;
      extremumIdx = cmpIdx;
      leftIdx = ((((cmpIdx - xstart) + 2) << 1) + xstart) - 1;
      changed = true;
    } else {
      exitg1 = true;
    }
  }

  if (changed && (leftIdx <= xend)) {
    extremum = x->data[extremumIdx];
    cmpIdx = x->data[leftIdx - 1];
    if (x->data[extremumIdx] < cmpIdx) {
      x->data[extremumIdx] = cmpIdx;
      x->data[leftIdx - 1] = extremum;
    }
  }
}

//
// Arguments    : emxArray_int32_T *x
//                int xstart
//                int xend
// Return Type  : void
//
void b_heapsort(emxArray_int32_T *x, int xstart, int xend)
{
  int n;
  int idx;
  int t;
  n = (xend - xstart) - 1;
  for (idx = n + 2; idx >= 1; idx--) {
    heapify(x, idx, xstart, xend);
  }

  for (idx = 0; idx <= n; idx++) {
    t = x->data[xend - 1];
    x->data[xend - 1] = x->data[xstart - 1];
    x->data[xstart - 1] = t;
    xend--;
    heapify(x, 1, xstart, xend);
  }
}

//
// File trailer for heapsort.cpp
//
// [EOF]
//
