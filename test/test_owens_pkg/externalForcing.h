//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: externalForcing.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//
#ifndef EXTERNALFORCING_H
#define EXTERNALFORCING_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void externalForcing(double time, const double aeroLoads_timeArray_data[],
  const int aeroLoads_timeArray_size[1], const emxArray_real_T
  *aeroLoads_ForceValHist, const emxArray_real_T *aeroLoads_ForceDof,
  emxArray_real_T *Fexternal, emxArray_real_T *Fdof);

#endif

//
// File trailer for externalForcing.h
//
// [EOF]
//
