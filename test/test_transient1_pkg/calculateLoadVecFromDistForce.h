/*
 * File: calculateLoadVecFromDistForce.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

#ifndef CALCULATELOADVECFROMDISTFORCE_H
#define CALCULATELOADVECFROMDISTFORCE_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "test_transient1_types.h"

/* Function Declarations */
extern void calculateLoadVecFromDistForce(const double input_xloc[2], const
  double input_sectionProps_twist[2], double input_sweepAngle, double
  input_coneAngle, double input_rollAngle, const double input_extDistF2Node[2],
  const double input_extDistF3Node[2], const double input_extDistF4Node[2],
  double output_Fe[12]);

#endif

/*
 * File trailer for calculateLoadVecFromDistForce.h
 *
 * [EOF]
 */
