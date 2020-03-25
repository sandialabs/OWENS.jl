/*
 * File: processAeroLoadsBLE.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

#ifndef PROCESSAEROLOADSBLE_H
#define PROCESSAEROLOADSBLE_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "test_transient1_types.h"

/* Function Declarations */
extern void processAeroLoadsBLE(const emxArray_char_T *aeroLoadsFile, const
  emxArray_char_T *OWENSfile, double aeroLoads_timeArray_data[], int
  aeroLoads_timeArray_size[1], emxArray_real_T *aeroLoads_ForceValHist,
  emxArray_real_T *aeroLoads_ForceDof);

#endif

/*
 * File trailer for processAeroLoadsBLE.h
 *
 * [EOF]
 */
