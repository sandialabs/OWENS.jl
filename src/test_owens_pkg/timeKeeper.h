//
// File: timeKeeper.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
//
#ifndef TIMEKEEPER_H
#define TIMEKEEPER_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void b_timeKeeper(double *outTime_tv_sec, double *outTime_tv_nsec);
extern void savedTime_not_empty_init();
extern void timeKeeper(double newTime_tv_sec, double newTime_tv_nsec);

#endif

//
// File trailer for timeKeeper.h
//
// [EOF]
//
