//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: sprintf.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:47:29
//

// Include Files
#include "sprintf.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <stdio.h>
#include <string.h>

// Function Definitions

//
// Arguments    : double varargin_1
//                emxArray_char_T *str
// Return Type  : void
//
void b_sprintf(double varargin_1, emxArray_char_T *str)
{
  int nbytes;
  int i;
  nbytes = snprintf(NULL, 0, "%1.12e", varargin_1);
  i = str->size[0] * str->size[1];
  str->size[0] = 1;
  str->size[1] = nbytes + 1;
  emxEnsureCapacity_char_T(str, i);
  snprintf(&str->data[0], (size_t)(nbytes + 1), "%1.12e", varargin_1);
  i = str->size[0] * str->size[1];
  if (1 > nbytes) {
    str->size[1] = 0;
  } else {
    str->size[1] = nbytes;
  }

  emxEnsureCapacity_char_T(str, i);
}

//
// File trailer for sprintf.cpp
//
// [EOF]
//
