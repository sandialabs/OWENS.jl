//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: find.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "find.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_boolean_T *x
//                emxArray_int32_T *i
// Return Type  : void
//
void b_eml_find(const emxArray_boolean_T *x, emxArray_int32_T *i)
{
  int nx;
  int idx;
  int ii;
  boolean_T exitg1;
  nx = x->size[0] * x->size[1];
  idx = 0;
  ii = i->size[0];
  i->size[0] = nx;
  emxEnsureCapacity_int32_T(i, ii);
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx - 1)) {
    if (x->data[ii]) {
      idx++;
      i->data[idx - 1] = ii + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }

  if (nx == 1) {
    if (idx == 0) {
      i->size[0] = 0;
    }
  } else {
    ii = i->size[0];
    if (1 > idx) {
      i->size[0] = 0;
    } else {
      i->size[0] = idx;
    }

    emxEnsureCapacity_int32_T(i, ii);
  }
}

//
// Arguments    : const emxArray_real_T *x
//                int i_data[]
//                int i_size[2]
// Return Type  : void
//
void c_eml_find(const emxArray_real_T *x, int i_data[], int i_size[2])
{
  int k;
  int idx;
  int ii;
  boolean_T exitg1;
  k = (1 <= x->size[1]);
  idx = 0;
  i_size[0] = 1;
  i_size[1] = k;
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= x->size[1] - 1)) {
    if (x->data[ii] != 0.0) {
      idx = 1;
      i_data[0] = ii + 1;
      exitg1 = true;
    } else {
      ii++;
    }
  }

  if (k == 1) {
    if (idx == 0) {
      i_size[0] = 1;
      i_size[1] = 0;
    }
  } else {
    i_size[1] = (1 <= idx);
  }
}

//
// Arguments    : const double x[8]
//                int i_data[]
//                int i_size[1]
// Return Type  : void
//
void d_eml_find(const double x[8], int i_data[], int i_size[1])
{
  int idx;
  int ii;
  boolean_T exitg1;
  idx = 0;
  i_size[0] = 1;
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii < 8)) {
    if (x[ii] != 0.0) {
      idx = 1;
      i_data[0] = ii + 1;
      exitg1 = true;
    } else {
      ii++;
    }
  }

  if (idx == 0) {
    i_size[0] = 0;
  }
}

//
// Arguments    : const emxArray_boolean_T *x
//                emxArray_int32_T *i
// Return Type  : void
//
void e_eml_find(const emxArray_boolean_T *x, emxArray_int32_T *i)
{
  int nx;
  int idx;
  int ii;
  boolean_T exitg1;
  nx = x->size[0];
  idx = 0;
  ii = i->size[0];
  i->size[0] = x->size[0];
  emxEnsureCapacity_int32_T(i, ii);
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx - 1)) {
    if (x->data[ii]) {
      idx++;
      i->data[idx - 1] = ii + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }

  if (x->size[0] == 1) {
    if (idx == 0) {
      i->size[0] = 0;
    }
  } else {
    ii = i->size[0];
    if (1 > idx) {
      i->size[0] = 0;
    } else {
      i->size[0] = idx;
    }

    emxEnsureCapacity_int32_T(i, ii);
  }
}

//
// Arguments    : int i_data[]
//                int i_size[2]
// Return Type  : void
//
void eml_find(int i_data[], int i_size[2])
{
  int idx;
  int ii;
  boolean_T exitg1;
  static const boolean_T b_bv[73] = { false, true, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, true, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false };

  idx = 0;
  i_size[0] = 1;
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii < 73)) {
    if (b_bv[ii]) {
      idx++;
      i_data[idx - 1] = ii + 1;
      if (idx >= 73) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }

  if (1 > idx) {
    i_size[1] = 0;
  } else {
    i_size[1] = idx;
  }
}

//
// File trailer for find.cpp
//
// [EOF]
//
