//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateReactionForceAtNode.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 03-Apr-2020 15:56:19
//
#ifndef CALCULATEREACTIONFORCEATNODE_H
#define CALCULATEREACTIONFORCEATNODE_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void calculateReactionForceAtNode(const char model_analysisType[3],
  double model_RayleighAlpha, double model_RayleighBeta, const emxArray_real_T
  *model_joint, const emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y,
  const emxArray_real_T *mesh_z, const emxArray_real_T *mesh_conn, const
  emxArray_struct_T *el_props, const emxArray_real_T *el_elLen, const
  emxArray_real_T *el_psi, const emxArray_real_T *el_theta, const
  emxArray_real_T *el_roll, const b_emxArray_struct_T *elStorage, const
  b_struct_T *timeInt, const i_struct_T dispData, const emxArray_real_T
  *displ_iter, double Omega, double OmegaDot, const double CN2H[9], double
  cummulativeForce[6]);

#endif

//
// File trailer for calculateReactionForceAtNode.h
//
// [EOF]
//
