//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateStrainForElements.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//
#ifndef CALCULATESTRAINFORELEMENTS_H
#define CALCULATESTRAINFORELEMENTS_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_transient1_types.h"

// Function Declarations
extern void calculateStrainForElements(double numEl, const emxArray_real_T *conn,
  const emxArray_struct_T *el_props, const emxArray_real_T *el_elLen, const
  emxArray_real_T *el_psi, const emxArray_real_T *el_theta, const
  emxArray_real_T *el_roll, const emxArray_real_T *displ, c_emxArray_struct_T
  *elStrain);

#endif

//
// File trailer for calculateStrainForElements.h
//
// [EOF]
//
