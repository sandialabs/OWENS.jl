/*
 * File: calculateShapeFunctions.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

#ifndef CALCULATESHAPEFUNCTIONS_H
#define CALCULATESHAPEFUNCTIONS_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "test_transient1_types.h"

/* Function Declarations */
extern void b_calculateShapeFunctions(const double x[2], double N_data[], int
  N_size[2], double p_N_x_data[], int p_N_x_size[1], double *Jac);
extern void calculateShapeFunctions(double xi, const double x[2], double N_data[],
  int N_size[2], double p_N_x_data[], int p_N_x_size[1], double *Jac);

#endif

/*
 * File trailer for calculateShapeFunctions.h
 *
 * [EOF]
 */
