//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateTimoshenkoElementNL.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//
#ifndef CALCULATETIMOSHENKOELEMENTNL_H
#define CALCULATETIMOSHENKOELEMENTNL_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void b_calculateTimoshenkoElementNL(const p_struct_T *input, const
  f_struct_T *elStorage, o_struct_T *output);
extern void c_calculateTimoshenkoElementNL(const double input_xloc[2], const
  double input_sectionProps_twist[2], const double input_sectionProps_rhoA[2],
  const double input_sectionProps_EA[2], const double input_sectionProps_zcm[2],
  const double input_sectionProps_ycm[2], double input_sweepAngle, double
  input_coneAngle, double input_rollAngle, const double input_concMass[8], const
  double input_disp_data[], const double input_x_data[], const double
  input_y_data[], const double input_z_data[], const double input_CN2H[9],
  double input_RayleighAlpha, double input_RayleighBeta, const f_struct_T
  *elStorage, o_struct_T *output);
extern void calculateTimoshenkoElementNL(const n_struct_T *input, const
  f_struct_T *elStorage, o_struct_T *output);

#endif

//
// File trailer for calculateTimoshenkoElementNL.h
//
// [EOF]
//
