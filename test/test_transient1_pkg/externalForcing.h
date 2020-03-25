/*
 * File: externalForcing.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

#ifndef EXTERNALFORCING_H
#define EXTERNALFORCING_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "test_transient1_types.h"

/* Function Declarations */
extern void externalForcing(double time, const double aeroLoads_timeArray_data[],
  const int aeroLoads_timeArray_size[1], const emxArray_real_T
  *aeroLoads_ForceValHist, const emxArray_real_T *aeroLoads_ForceDof,
  emxArray_real_T *Fexternal, emxArray_real_T *Fdof);

#endif

/*
 * File trailer for externalForcing.h
 *
 * [EOF]
 */
