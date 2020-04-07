//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: eye.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "eye.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// Arguments    : double b_I[9]
// Return Type  : void
//
void eye(double b_I[9])
{
  std::memset(&b_I[0], 0, 9U * sizeof(double));
  b_I[0] = 1.0;
  b_I[4] = 1.0;
  b_I[8] = 1.0;
}

//
// File trailer for eye.cpp
//
// [EOF]
//
