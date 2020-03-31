//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: mapCactusLoadsFile.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//
#ifndef MAPCACTUSLOADSFILE_H
#define MAPCACTUSLOADSFILE_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_transient1_types.h"

// Function Declarations
extern void mapCactusLoadsFile(const emxArray_char_T *geomFn, const
  emxArray_char_T *loadsFn, const emxArray_char_T *bldFn, const emxArray_char_T *
  elFn, const emxArray_char_T *ortFn, const emxArray_char_T *meshFn, double
  time_data[], int time_size[1], emxArray_real_T *ForceValHist, emxArray_real_T *
  ForceDof);

#endif

//
// File trailer for mapCactusLoadsFile.h
//
// [EOF]
//
