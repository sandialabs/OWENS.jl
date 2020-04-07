//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: writeOutput.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 15:21:39
//
#ifndef WRITEOUTPUT_H
#define WRITEOUTPUT_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void writeOutput(const emxArray_real_T *freq, const emxArray_real_T *damp,
  const emxArray_real_T *phase1, const emxArray_real_T *phase2, const
  emxArray_real_T *imagComponentSign, double fid, emxArray_real_T *freqSorted,
  emxArray_real_T *dampSorted, emxArray_real_T *imagCompSignSorted);

#endif

//
// File trailer for writeOutput.h
//
// [EOF]
//
