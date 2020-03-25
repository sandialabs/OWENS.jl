/*
 * File: elementPostProcess.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

#ifndef ELEMENTPOSTPROCESS_H
#define ELEMENTPOSTPROCESS_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "test_transient1_types.h"

/* Function Declarations */
extern void elementPostProcess(double elementNumber, const char
  model_analysisType[3], double model_RayleighAlpha, double model_RayleighBeta,
  const emxArray_real_T *model_joint, const emxArray_real_T *mesh_x, const
  emxArray_real_T *mesh_y, const emxArray_real_T *mesh_z, const emxArray_real_T *
  mesh_conn, const emxArray_struct_T *el_props, const emxArray_real_T *el_elLen,
  const emxArray_real_T *el_psi, const emxArray_real_T *el_theta, const
  emxArray_real_T *el_roll, const c_emxArray_struct_T *elStorage, const struct_T
  *timeInt, const g_struct_T dispData, const emxArray_real_T *displ_iter, const
  double rbData[9], double Omega, double OmegaDot, const double CN2H[9], double
  Fpp[12]);

#endif

/*
 * File trailer for elementPostProcess.h
 *
 * [EOF]
 */
