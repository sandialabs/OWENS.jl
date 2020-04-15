//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: staticAnalysis.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//
#ifndef STATICANALYSIS_H
#define STATICANALYSIS_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void b_applyConstraints(emxArray_real_T *Kg, const emxArray_real_T
  *transMatrix);
extern void staticAnalysis(const char model_nlParams_iterationType[2], double
  model_RayleighAlpha, double model_RayleighBeta, double model_BC_numpBC, const
  emxArray_real_T *model_BC_pBC, const emxArray_real_T *model_joint, const
  emxArray_real_T *model_jointTransform, double mesh_numEl, const
  emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y, const emxArray_real_T *
  mesh_z, const emxArray_real_T *mesh_conn, const i_struct_T el, emxArray_real_T
  *displ, const c_emxArray_struct_T *elStorage, b_emxArray_struct_T *elStrain,
  boolean_T *staticAnalysisSuccessful);

#endif

//
// File trailer for staticAnalysis.h
//
// [EOF]
//
