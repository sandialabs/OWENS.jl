//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: timeKeeper.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "timeKeeper.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Type Definitions
struct struct_T
{
  double tv_sec;
  double tv_nsec;
};

// Variable Definitions
static struct_T savedTime;
static boolean_T savedTime_not_empty;

// Function Definitions

//
// Arguments    : double *outTime_tv_sec
//                double *outTime_tv_nsec
// Return Type  : void
//
void b_timeKeeper(double *outTime_tv_sec, double *outTime_tv_nsec)
{
  *outTime_tv_sec = savedTime.tv_sec;
  *outTime_tv_nsec = savedTime.tv_nsec;
}

//
// Arguments    : void
// Return Type  : void
//
void savedTime_not_empty_init()
{
  savedTime_not_empty = false;
}

//
// Arguments    : double newTime_tv_sec
//                double newTime_tv_nsec
// Return Type  : void
//
void timeKeeper(double newTime_tv_sec, double newTime_tv_nsec)
{
  struct timespec b_timespec;
  if (!savedTime_not_empty) {
    clock_gettime(CLOCK_MONOTONIC, &b_timespec);
    savedTime_not_empty = true;
  }

  savedTime.tv_sec = newTime_tv_sec;
  savedTime.tv_nsec = newTime_tv_nsec;
}

//
// File trailer for timeKeeper.cpp
//
// [EOF]
//
