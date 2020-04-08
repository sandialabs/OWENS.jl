//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateLoadVecFromDistForce.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:47:29
//
#ifndef CALCULATELOADVECFROMDISTFORCE_H
#define CALCULATELOADVECFROMDISTFORCE_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void calculateLoadVecFromDistForce(const double input_xloc[2], const
  double input_sectionProps_twist[2], double input_sweepAngle, double
  input_coneAngle, double input_rollAngle, const double input_extDistF2Node[2],
  const double input_extDistF3Node[2], const double input_extDistF4Node[2],
  double output_Fe[12]);

#endif

//
// File trailer for calculateLoadVecFromDistForce.h
//
// [EOF]
//
