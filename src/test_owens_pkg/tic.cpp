//
// File: tic.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
//

// Include Files
#include "tic.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
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
  timeKeeper((double)b_timespec.tv_sec, (double)b_timespec.tv_nsec);
}

//
// File trailer for tic.cpp
//
// [EOF]
//
