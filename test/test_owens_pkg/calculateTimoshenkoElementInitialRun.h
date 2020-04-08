//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateTimoshenkoElementInitialRun.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:47:29
//
#ifndef CALCULATETIMOSHENKOELEMENTINITIALRUN_H
#define CALCULATETIMOSHENKOELEMENTINITIALRUN_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void c_calculateTimoshenkoElementIni(const double input_xloc[2], const
  double input_sectionProps_twist[2], const double input_sectionProps_rhoA[2],
  const double input_sectionProps_EIyy[2], const double input_sectionProps_EIzz
  [2], const double input_sectionProps_GJ[2], const double
  input_sectionProps_EA[2], const double input_sectionProps_rhoIyy[2], const
  double input_sectionProps_rhoIzz[2], const double input_sectionProps_rhoJ[2],
  const double input_sectionProps_zcm[2], const double input_sectionProps_ycm[2],
  double input_sweepAngle, double input_coneAngle, double input_rollAngle, const
  double input_x[2], const double input_y[2], const double input_z[2], boolean_T
  input_concMassFlag, const double input_concMass[8], f_struct_T *elStorage);

#endif

//
// File trailer for calculateTimoshenkoElementInitialRun.h
//
// [EOF]
//
