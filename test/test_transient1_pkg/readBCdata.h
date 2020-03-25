/*
 * File: readBCdata.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

#ifndef READBCDATA_H
#define READBCDATA_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "test_transient1_types.h"

/* Function Declarations */
extern void readBCdata(const emxArray_char_T *bcfilename, double numNodes,
  double *BC_numpBC, emxArray_real_T *BC_pBC, double *BC_numsBC, double
  *BC_nummBC, emxArray_real_T *BC_isConstrained);

#endif

/*
 * File trailer for readBCdata.h
 *
 * [EOF]
 */
