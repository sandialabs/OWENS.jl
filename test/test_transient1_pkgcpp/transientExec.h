//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: transientExec.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//
#ifndef TRANSIENTEXEC_H
#define TRANSIENTEXEC_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_transient1_types.h"

// Function Declarations
extern void transientExec(const char model_analysisType[3], const double
  model_tocp[2], const emxArray_char_T *model_aeroloadfile, const
  emxArray_char_T *model_owensfile, double model_RayleighAlpha, double
  model_RayleighBeta, double model_BC_numpBC, const emxArray_real_T
  *model_BC_pBC, const emxArray_real_T *model_joint, const emxArray_real_T
  *model_jointTransform, const e_struct_T mesh, const f_struct_T el);

#endif

//
// File trailer for transientExec.h
//
// [EOF]
//
