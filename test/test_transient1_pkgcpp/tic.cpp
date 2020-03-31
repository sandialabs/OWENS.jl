//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: tic.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "tic.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "timeKeeper.h"
#include <string.h>

// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void tic()
{
  struct timespec b_timespec;
  clock_gettime(CLOCK_MONOTONIC, &b_timespec);
  timeKeeper();
}

//
// File trailer for tic.cpp
//
// [EOF]
//
