//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: initialElementCalculations.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//
#ifndef INITIALELEMENTCALCULATIONS_H
#define INITIALELEMENTCALCULATIONS_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void initialElementCalculations(const emxArray_real_T *model_joint, const
  emxArray_struct_T *el_props, const emxArray_real_T *el_elLen, const
  emxArray_real_T *el_psi, const emxArray_real_T *el_theta, const
  emxArray_real_T *el_roll, double mesh_numEl, const emxArray_real_T *mesh_x,
  const emxArray_real_T *mesh_y, const emxArray_real_T *mesh_z, const
  emxArray_real_T *mesh_conn, c_emxArray_struct_T *elStorage);

#endif

//
// File trailer for initialElementCalculations.h
//
// [EOF]
//
