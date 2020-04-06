//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: extractFreqDamp.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//
#ifndef EXTRACTFREQDAMP_H
#define EXTRACTFREQDAMP_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void extractFreqDamp(const creal_T val, const emxArray_creal_T *vec,
  const emxArray_real_T *jointTransform, const emxArray_real_T *reducedDOFList,
  double BC_numpBC, const emxArray_real_T *BC_pBC, double *freq, double *damp,
  emxArray_real_T *phase1, emxArray_real_T *phase2, emxArray_creal_T
  *sortedModes);

#endif

//
// File trailer for extractFreqDamp.h
//
// [EOF]
//
