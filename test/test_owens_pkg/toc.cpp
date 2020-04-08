//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: toc.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:47:29
//

// Include Files
#include "toc.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "timeKeeper.h"
#include <string.h>

// Function Definitions

//
// Arguments    : void
// Return Type  : double
//
double toc()
{
  double tstart_tv_sec;
  double tstart_tv_nsec;
  struct timespec b_timespec;
  b_timeKeeper(&tstart_tv_sec, &tstart_tv_nsec);
  clock_gettime(CLOCK_MONOTONIC, &b_timespec);
  return ((double)b_timespec.tv_sec - tstart_tv_sec) + ((double)
    b_timespec.tv_nsec - tstart_tv_nsec) / 1.0E+9;
}

//
// File trailer for toc.cpp
//
// [EOF]
//
