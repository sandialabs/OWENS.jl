//
// File: readNLParamsFile.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
//
#ifndef READNLPARAMSFILE_H
#define READNLPARAMSFILE_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void readNLParamsFile(char nlParams_iterationType[2], boolean_T
  *c_nlParams_adaptiveLoadStepping, double *nlParams_tolerance, double
  *nlParams_maxIterations, double *nlParams_maxNumLoadSteps, double
  *nlParams_minLoadStepDelta, double *nlParams_minLoadStep, double
  *nlParams_prescribedLoadStep);

#endif

//
// File trailer for readNLParamsFile.h
//
// [EOF]
//
