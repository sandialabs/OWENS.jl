//
// File: structuralDynamicsTransient.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
//
#ifndef STRUCTURALDYNAMICSTRANSIENT_H
#define STRUCTURALDYNAMICSTRANSIENT_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void structuralDynamicsTransient(const char model_analysisType[3], double
  model_RayleighAlpha, double model_RayleighBeta, double model_BC_numpBC, const
  emxArray_real_T *model_BC_pBC, const emxArray_real_T *model_joint, const
  emxArray_real_T *model_jointTransform, double mesh_numEl, const
  emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y, const emxArray_real_T *
  mesh_z, const emxArray_real_T *mesh_conn, const i_struct_T el, const
  j_struct_T dispData, double Omega, double OmegaDot, const c_emxArray_struct_T *
  elStorage, const emxArray_real_T *Fexternal, const emxArray_real_T *Fdof,
  const double CN2H[9], k_struct_T *dispOut, double FReaction_sp1[6]);

#endif

//
// File trailer for structuralDynamicsTransient.h
//
// [EOF]
//
