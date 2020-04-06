//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readBCdata.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//
#ifndef READBCDATA_H
#define READBCDATA_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void readBCdata(const emxArray_char_T *bcfilename, double numNodes,
  double *BC_numpBC, emxArray_real_T *BC_pBC, double *BC_numsBC, double
  *BC_nummBC, emxArray_real_T *BC_isConstrained);

#endif

//
// File trailer for readBCdata.h
//
// [EOF]
//
