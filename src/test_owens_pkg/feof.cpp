//
// File: feof.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
//

// Include Files
#include "feof.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// Arguments    : double fileID
// Return Type  : double
//
double b_feof(double fileID)
{
  double st;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T a;
  int b_st;
  b_NULL = NULL;
  getfilestar(fileID, &filestar, &a);
  if (filestar == b_NULL) {
    st = 0.0;
  } else {
    b_st = feof(filestar);
    st = ((int)b_st != 0);
  }

  return st;
}

//
// File trailer for feof.cpp
//
// [EOF]
//
