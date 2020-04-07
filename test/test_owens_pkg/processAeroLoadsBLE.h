//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: processAeroLoadsBLE.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//
#ifndef PROCESSAEROLOADSBLE_H
#define PROCESSAEROLOADSBLE_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void processAeroLoadsBLE(const emxArray_char_T *aeroLoadsFile, const
  emxArray_char_T *OWENSfile, double aeroLoads_timeArray_data[], int
  aeroLoads_timeArray_size[1], emxArray_real_T *aeroLoads_ForceValHist,
  emxArray_real_T *aeroLoads_ForceDof);

#endif

//
// File trailer for processAeroLoadsBLE.h
//
// [EOF]
//
