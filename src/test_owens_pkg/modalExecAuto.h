//
// File: modalExecAuto.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//
#ifndef MODALEXECAUTO_H
#define MODALEXECAUTO_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void modalExecAuto(const char model_analysisType[2], const char
  model_nlParams_iterationType[2], double model_RayleighAlpha, double
  model_RayleighBeta, double model_BC_numpBC, const emxArray_real_T
  *model_BC_pBC, const emxArray_real_T *model_BC_map, const emxArray_real_T
  *model_joint, const char model_outFilename_data[], const int
  model_outFilename_size[2], const emxArray_real_T *model_jointTransform, const
  emxArray_real_T *model_reducedDOFList, const h_struct_T mesh, const i_struct_T
  el, const emxArray_real_T *displ, emxArray_real_T *freq, emxArray_real_T *damp);

#endif

//
// File trailer for modalExecAuto.h
//
// [EOF]
//
